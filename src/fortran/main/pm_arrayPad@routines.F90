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
!>  This file contains procedure implementations of [pm_arrayPad](@ref pm_arrayPad).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_arrayPad) routines ! LCOV_EXCL_LINE

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

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPadded_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Asis_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getPaddedAsisSB_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPaddedAsisSB_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPaddedAsisSB_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPaddedAsisSB_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPaddedAsisSB_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getPaddedAsisSB_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPaddedAsisSB_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPaddedAsisSB_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPaddedAsisSB_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPaddedAsisSB_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getPaddedAsisSB_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getPaddedAsisSB_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getPaddedAsisSB_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getPaddedAsisSB_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getPaddedAsisSB_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getPaddedAsisSB_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getPaddedAsisSB_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getPaddedAsisSB_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getPaddedAsisSB_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getPaddedAsisSB_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPaddedAsisSB_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPaddedAsisSB_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPaddedAsisSB_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPaddedAsisSB_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPaddedAsisSB_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPaddedAsisSB_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPaddedAsisSB_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPaddedAsisSB_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPaddedAsisSB_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPaddedAsisSB_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Asis_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Marg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getPaddedMargSB_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPaddedMargSB_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPaddedMargSB_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPaddedMargSB_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPaddedMargSB_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getPaddedMargSB_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPaddedMargSB_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPaddedMargSB_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPaddedMargSB_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPaddedMargSB_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getPaddedMargSB_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getPaddedMargSB_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getPaddedMargSB_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getPaddedMargSB_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getPaddedMargSB_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getPaddedMargSB_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getPaddedMargSB_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getPaddedMargSB_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getPaddedMargSB_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getPaddedMargSB_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPaddedMargSB_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPaddedMargSB_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPaddedMargSB_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPaddedMargSB_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPaddedMargSB_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPaddedMargSB_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPaddedMargSB_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPaddedMargSB_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPaddedMargSB_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPaddedMargSB_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Marg_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPadded_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPadded_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Asis_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setPaddedAsisSB_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setPaddedAsisSB_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setPaddedAsisSB_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setPaddedAsisSB_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setPaddedAsisSB_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setPaddedAsisSB_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setPaddedAsisSB_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setPaddedAsisSB_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setPaddedAsisSB_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setPaddedAsisSB_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setPaddedAsisSB_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setPaddedAsisSB_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setPaddedAsisSB_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setPaddedAsisSB_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setPaddedAsisSB_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setPaddedAsisSB_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setPaddedAsisSB_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setPaddedAsisSB_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setPaddedAsisSB_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setPaddedAsisSB_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setPaddedAsisSB_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPaddedAsisSB_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPaddedAsisSB_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPaddedAsisSB_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPaddedAsisSB_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPaddedAsisSB_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPaddedAsisSB_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPaddedAsisSB_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPaddedAsisSB_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPaddedAsisSB_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Asis_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Marg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setPaddedMargSB_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setPaddedMargSB_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setPaddedMargSB_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setPaddedMargSB_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setPaddedMargSB_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setPaddedMargSB_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setPaddedMargSB_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setPaddedMargSB_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setPaddedMargSB_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setPaddedMargSB_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setPaddedMargSB_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setPaddedMargSB_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setPaddedMargSB_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setPaddedMargSB_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setPaddedMargSB_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setPaddedMargSB_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setPaddedMargSB_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setPaddedMargSB_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setPaddedMargSB_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setPaddedMargSB_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setPaddedMargSB_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPaddedMargSB_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPaddedMargSB_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPaddedMargSB_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPaddedMargSB_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPaddedMargSB_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPaddedMargSB_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPaddedMargSB_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPaddedMargSB_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPaddedMargSB_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Marg_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPadded_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPadded_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Asis_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getPaddedAsisSL_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPaddedAsisSL_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPaddedAsisSL_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPaddedAsisSL_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPaddedAsisSL_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getPaddedAsisSL_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPaddedAsisSL_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPaddedAsisSL_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPaddedAsisSL_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPaddedAsisSL_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getPaddedAsisSL_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getPaddedAsisSL_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getPaddedAsisSL_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getPaddedAsisSL_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getPaddedAsisSL_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getPaddedAsisSL_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getPaddedAsisSL_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getPaddedAsisSL_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getPaddedAsisSL_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getPaddedAsisSL_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPaddedAsisSL_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPaddedAsisSL_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPaddedAsisSL_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPaddedAsisSL_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPaddedAsisSL_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPaddedAsisSL_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPaddedAsisSL_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPaddedAsisSL_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPaddedAsisSL_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPaddedAsisSL_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Asis_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Marg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getPaddedMargSL_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPaddedMargSL_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPaddedMargSL_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPaddedMargSL_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPaddedMargSL_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getPaddedMargSL_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPaddedMargSL_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPaddedMargSL_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPaddedMargSL_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPaddedMargSL_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getPaddedMargSL_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getPaddedMargSL_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getPaddedMargSL_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getPaddedMargSL_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getPaddedMargSL_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getPaddedMargSL_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getPaddedMargSL_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getPaddedMargSL_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getPaddedMargSL_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getPaddedMargSL_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPaddedMargSL_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPaddedMargSL_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPaddedMargSL_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPaddedMargSL_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPaddedMargSL_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPaddedMargSL_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPaddedMargSL_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPaddedMargSL_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPaddedMargSL_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPaddedMargSL_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Marg_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPadded_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPadded_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Asis_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setPaddedAsisSL_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setPaddedAsisSL_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setPaddedAsisSL_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setPaddedAsisSL_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setPaddedAsisSL_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setPaddedAsisSL_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setPaddedAsisSL_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setPaddedAsisSL_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setPaddedAsisSL_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setPaddedAsisSL_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setPaddedAsisSL_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setPaddedAsisSL_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setPaddedAsisSL_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setPaddedAsisSL_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setPaddedAsisSL_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setPaddedAsisSL_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setPaddedAsisSL_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setPaddedAsisSL_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setPaddedAsisSL_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setPaddedAsisSL_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setPaddedAsisSL_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPaddedAsisSL_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPaddedAsisSL_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPaddedAsisSL_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPaddedAsisSL_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPaddedAsisSL_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPaddedAsisSL_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPaddedAsisSL_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPaddedAsisSL_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPaddedAsisSL_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Asis_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Marg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setPaddedMargSL_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setPaddedMargSL_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setPaddedMargSL_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setPaddedMargSL_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setPaddedMargSL_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setPaddedMargSL_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setPaddedMargSL_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setPaddedMargSL_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setPaddedMargSL_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setPaddedMargSL_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setPaddedMargSL_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setPaddedMargSL_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setPaddedMargSL_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setPaddedMargSL_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setPaddedMargSL_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setPaddedMargSL_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setPaddedMargSL_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setPaddedMargSL_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setPaddedMargSL_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setPaddedMargSL_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setPaddedMargSL_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPaddedMargSL_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPaddedMargSL_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPaddedMargSL_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPaddedMargSL_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPaddedMargSL_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPaddedMargSL_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPaddedMargSL_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPaddedMargSL_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPaddedMargSL_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Marg_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPadded_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPadded_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Asis_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getPaddedAsisSR_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPaddedAsisSR_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPaddedAsisSR_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPaddedAsisSR_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPaddedAsisSR_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getPaddedAsisSR_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPaddedAsisSR_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPaddedAsisSR_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPaddedAsisSR_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPaddedAsisSR_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getPaddedAsisSR_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getPaddedAsisSR_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getPaddedAsisSR_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getPaddedAsisSR_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getPaddedAsisSR_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getPaddedAsisSR_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getPaddedAsisSR_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getPaddedAsisSR_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getPaddedAsisSR_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getPaddedAsisSR_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPaddedAsisSR_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPaddedAsisSR_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPaddedAsisSR_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPaddedAsisSR_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPaddedAsisSR_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPaddedAsisSR_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPaddedAsisSR_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPaddedAsisSR_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPaddedAsisSR_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPaddedAsisSR_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Asis_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Marg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getPaddedMargSR_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPaddedMargSR_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPaddedMargSR_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPaddedMargSR_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPaddedMargSR_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getPaddedMargSR_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPaddedMargSR_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPaddedMargSR_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPaddedMargSR_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPaddedMargSR_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getPaddedMargSR_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getPaddedMargSR_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getPaddedMargSR_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getPaddedMargSR_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getPaddedMargSR_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getPaddedMargSR_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getPaddedMargSR_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getPaddedMargSR_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getPaddedMargSR_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getPaddedMargSR_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPaddedMargSR_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPaddedMargSR_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPaddedMargSR_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPaddedMargSR_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPaddedMargSR_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPaddedMargSR_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPaddedMargSR_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPaddedMargSR_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPaddedMargSR_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPaddedMargSR_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Marg_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPadded_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPadded_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Asis_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setPaddedAsisSR_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setPaddedAsisSR_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setPaddedAsisSR_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setPaddedAsisSR_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setPaddedAsisSR_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setPaddedAsisSR_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setPaddedAsisSR_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setPaddedAsisSR_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setPaddedAsisSR_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setPaddedAsisSR_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setPaddedAsisSR_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setPaddedAsisSR_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setPaddedAsisSR_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setPaddedAsisSR_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setPaddedAsisSR_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setPaddedAsisSR_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setPaddedAsisSR_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setPaddedAsisSR_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setPaddedAsisSR_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setPaddedAsisSR_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setPaddedAsisSR_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPaddedAsisSR_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPaddedAsisSR_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPaddedAsisSR_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPaddedAsisSR_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPaddedAsisSR_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPaddedAsisSR_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPaddedAsisSR_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPaddedAsisSR_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPaddedAsisSR_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Asis_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Marg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setPaddedMargSR_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setPaddedMargSR_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setPaddedMargSR_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setPaddedMargSR_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setPaddedMargSR_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setPaddedMargSR_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setPaddedMargSR_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setPaddedMargSR_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setPaddedMargSR_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setPaddedMargSR_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setPaddedMargSR_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setPaddedMargSR_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setPaddedMargSR_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setPaddedMargSR_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setPaddedMargSR_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setPaddedMargSR_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setPaddedMargSR_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setPaddedMargSR_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setPaddedMargSR_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setPaddedMargSR_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setPaddedMargSR_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPaddedMargSR_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPaddedMargSR_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPaddedMargSR_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPaddedMargSR_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPaddedMargSR_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPaddedMargSR_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPaddedMargSR_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPaddedMargSR_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPaddedMargSR_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Marg_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPadded_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines