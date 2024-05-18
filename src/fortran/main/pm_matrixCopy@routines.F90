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
!>  This file contains procedure implementations of [pm_matrixCopy](@ref pm_matrixCopy).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_matrixCopy) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_mathSqrt, only: getSqrt
    use pm_matrixIndex, only: getMatIndex
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatCopy_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RDP_RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXX_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXX_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXX_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXX_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXX_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXX_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXX_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLX_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLX_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLX_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLX_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLX_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLX_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLX_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XXD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_UXD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_XLD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULX_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULX_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULX_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULX_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULX_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RDP_ULX_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULX_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_RDP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define RDP_LFP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XXD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_LFP_UXD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_LFP_XLD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_LFP_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED || CHECK_ENABLED
#define LFP_RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XXD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_LFP_RDP_UXD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_LFP_RDP_XLD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LFP_RDP_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define RDP_RFP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RFP_UXD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RDP_RFP_XLD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_RFP_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define RFP_RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RFP_RDP_UXD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatCopy_RFP_RDP_XLD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RFP_RDP_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMatCopy_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatCopy_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RDP_RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXX_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXX_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXX_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXX_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXX_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXX_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXX_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLX_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLX_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLX_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLX_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLX_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLX_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLX_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XXD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_UXD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_XLD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULX_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULX_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULX_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULX_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULX_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RDP_ULX_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULX_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_RDP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define RDP_LFP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XXD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_LFP_UXD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_LFP_XLD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_LFP_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED || CHECK_ENABLED
#define LFP_RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XXD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_LFP_RDP_UXD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define TSO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_TSO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef TSO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define THO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_LFP_RDP_XLD_THO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef THO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LFP_RDP_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define RDP_RFP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RFP_UXD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RDP_RFP_XLD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_RFP_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define RFP_RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RFP_RDP_UXD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1
#define AIO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatCopy_RFP_RDP_XLD_AIO_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixCopy@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED
#undef AIO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RFP_RDP_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatCopy_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines
