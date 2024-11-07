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
!>  This file contains procedure implementations of [pm_arraySort](@ref pm_arraySort).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_arraySort) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

#if PDT_ENABLED
    use pm_container, only: css_pdt
#endif
    use pm_container, only: css_type
    use pm_arrayMerge, only: setMerged

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isAscendingAll_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isAscendingAllDefCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isAscendingAllDefCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isAscendingAllDefCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isAscendingAllDefCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isAscendingAllDefCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure isAscendingAllDefCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isAscendingAllDefCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isAscendingAllDefCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isAscendingAllDefCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isAscendingAllDefCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure isAscendingAllDefCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure isAscendingAllDefCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure isAscendingAllDefCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure isAscendingAllDefCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure isAscendingAllDefCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure isAscendingAllDefCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure isAscendingAllDefCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure isAscendingAllDefCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure isAscendingAllDefCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure isAscendingAllDefCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isAscendingAllDefCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isAscendingAllDefCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isAscendingAllDefCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isAscendingAllDefCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isAscendingAllDefCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isAscendingAllDefCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isAscendingAllDefCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isAscendingAllDefCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isAscendingAllDefCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isAscendingAllDefCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure isAscendingAllDefCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isAscendingAllDefCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isAscendingAllDefCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isAscendingAllDefCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isAscendingAllDefCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure isAscendingAllDefCom_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isAscendingAll_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isDescendingAll_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isDescendingAllDefCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isDescendingAllDefCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isDescendingAllDefCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isDescendingAllDefCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isDescendingAllDefCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure isDescendingAllDefCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isDescendingAllDefCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isDescendingAllDefCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isDescendingAllDefCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isDescendingAllDefCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure isDescendingAllDefCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure isDescendingAllDefCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure isDescendingAllDefCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure isDescendingAllDefCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure isDescendingAllDefCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure isDescendingAllDefCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure isDescendingAllDefCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure isDescendingAllDefCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure isDescendingAllDefCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure isDescendingAllDefCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isDescendingAllDefCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isDescendingAllDefCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isDescendingAllDefCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isDescendingAllDefCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isDescendingAllDefCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isDescendingAllDefCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isDescendingAllDefCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isDescendingAllDefCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isDescendingAllDefCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isDescendingAllDefCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure isDescendingAllDefCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isDescendingAllDefCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isDescendingAllDefCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isDescendingAllDefCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isDescendingAllDefCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure isDescendingAllDefCom_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isDescendingAll_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isAscending_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isAscendingDefCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isAscendingDefCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isAscendingDefCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isAscendingDefCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isAscendingDefCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure isAscendingDefCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isAscendingDefCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isAscendingDefCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isAscendingDefCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isAscendingDefCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure isAscendingDefCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure isAscendingDefCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure isAscendingDefCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure isAscendingDefCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure isAscendingDefCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure isAscendingDefCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure isAscendingDefCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure isAscendingDefCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure isAscendingDefCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure isAscendingDefCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isAscendingDefCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isAscendingDefCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isAscendingDefCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isAscendingDefCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isAscendingDefCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isAscendingDefCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isAscendingDefCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isAscendingDefCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isAscendingDefCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isAscendingDefCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure isAscendingDefCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isAscendingDefCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isAscendingDefCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isAscendingDefCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isAscendingDefCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure isAscendingDefCom_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isAscending_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isDescending_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isDescendingDefCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isDescendingDefCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isDescendingDefCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isDescendingDefCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isDescendingDefCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure isDescendingDefCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isDescendingDefCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isDescendingDefCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isDescendingDefCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isDescendingDefCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure isDescendingDefCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure isDescendingDefCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure isDescendingDefCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure isDescendingDefCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure isDescendingDefCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure isDescendingDefCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure isDescendingDefCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure isDescendingDefCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure isDescendingDefCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure isDescendingDefCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isDescendingDefCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isDescendingDefCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isDescendingDefCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isDescendingDefCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isDescendingDefCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isDescendingDefCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isDescendingDefCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isDescendingDefCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isDescendingDefCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isDescendingDefCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure isDescendingDefCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isDescendingDefCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isDescendingDefCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isDescendingDefCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isDescendingDefCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure isDescendingDefCom_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isDescending_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isSortedDefCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isSortedDefCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isSortedDefCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isSortedDefCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isSortedDefCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure isSortedDefCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isSortedDefCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isSortedDefCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isSortedDefCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isSortedDefCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure isSortedDefCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure isSortedDefCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure isSortedDefCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure isSortedDefCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure isSortedDefCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure isSortedDefCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure isSortedDefCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure isSortedDefCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure isSortedDefCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure isSortedDefCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isSortedDefCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isSortedDefCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isSortedDefCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isSortedDefCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isSortedDefCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isSortedDefCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isSortedDefCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isSortedDefCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isSortedDefCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isSortedDefCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure isSortedDefCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isSortedDefCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isSortedDefCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isSortedDefCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isSortedDefCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure isSortedDefCom_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isSortedCusCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isSortedCusCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isSortedCusCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isSortedCusCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isSortedCusCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure isSortedCusCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isSortedCusCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isSortedCusCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isSortedCusCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isSortedCusCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure isSortedCusCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure isSortedCusCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure isSortedCusCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure isSortedCusCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure isSortedCusCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure isSortedCusCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure isSortedCusCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure isSortedCusCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure isSortedCusCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure isSortedCusCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isSortedCusCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isSortedCusCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isSortedCusCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isSortedCusCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isSortedCusCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isSortedCusCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isSortedCusCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isSortedCusCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isSortedCusCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isSortedCusCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure isSortedCusCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isSortedCusCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isSortedCusCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isSortedCusCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isSortedCusCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure isSortedCusCom_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ind_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getSortedIndCusComDef_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getSortedIndCusComDef_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getSortedIndCusComDef_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getSortedIndCusComDef_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getSortedIndCusComDef_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure getSortedIndCusComDef_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getSortedIndCusComDef_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getSortedIndCusComDef_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getSortedIndCusComDef_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getSortedIndCusComDef_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getSortedIndCusComDef_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getSortedIndCusComDef_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getSortedIndCusComDef_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getSortedIndCusComDef_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getSortedIndCusComDef_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getSortedIndCusComDef_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getSortedIndCusComDef_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getSortedIndCusComDef_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getSortedIndCusComDef_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getSortedIndCusComDef_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getSortedIndCusComDef_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getSortedIndCusComDef_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getSortedIndCusComDef_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getSortedIndCusComDef_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getSortedIndCusComDef_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getSortedIndCusComDef_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getSortedIndCusComDef_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getSortedIndCusComDef_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getSortedIndCusComDef_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getSortedIndCusComDef_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getSortedIndCusComDef_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getSortedIndCusComDef_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getSortedIndCusComDef_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getSortedIndCusComDef_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getSortedIndCusComDef_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getSortedIndCusComDef_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ind_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getSortedArrCusComDef_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getSortedArrCusComDef_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getSortedArrCusComDef_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getSortedArrCusComDef_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getSortedArrCusComDef_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure getSortedArrCusComDef_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getSortedArrCusComDef_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getSortedArrCusComDef_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getSortedArrCusComDef_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getSortedArrCusComDef_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getSortedArrCusComDef_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getSortedArrCusComDef_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getSortedArrCusComDef_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getSortedArrCusComDef_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getSortedArrCusComDef_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getSortedArrCusComDef_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getSortedArrCusComDef_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getSortedArrCusComDef_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getSortedArrCusComDef_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getSortedArrCusComDef_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getSortedArrCusComDef_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getSortedArrCusComDef_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getSortedArrCusComDef_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getSortedArrCusComDef_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getSortedArrCusComDef_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getSortedArrCusComDef_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getSortedArrCusComDef_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getSortedArrCusComDef_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getSortedArrCusComDef_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getSortedArrCusComDef_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getSortedArrCusComDef_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getSortedArrCusComDef_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getSortedArrCusComDef_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getSortedArrCusComDef_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getSortedArrCusComDef_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getSortedArrCusComDef_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ind_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getSortedIndDefComDef_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getSortedIndDefComDef_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getSortedIndDefComDef_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getSortedIndDefComDef_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getSortedIndDefComDef_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure getSortedIndDefComDef_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getSortedIndDefComDef_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getSortedIndDefComDef_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getSortedIndDefComDef_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getSortedIndDefComDef_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getSortedIndDefComDef_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getSortedIndDefComDef_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getSortedIndDefComDef_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getSortedIndDefComDef_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getSortedIndDefComDef_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getSortedIndDefComDef_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getSortedIndDefComDef_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getSortedIndDefComDef_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getSortedIndDefComDef_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getSortedIndDefComDef_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getSortedIndDefComDef_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getSortedIndDefComDef_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getSortedIndDefComDef_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getSortedIndDefComDef_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getSortedIndDefComDef_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getSortedIndDefComDef_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getSortedIndDefComDef_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getSortedIndDefComDef_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getSortedIndDefComDef_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getSortedIndDefComDef_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getSortedIndDefComDef_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getSortedIndDefComDef_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getSortedIndDefComDef_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getSortedIndDefComDef_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getSortedIndDefComDef_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getSortedIndDefComDef_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ind_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getSortedArrDefComDef_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getSortedArrDefComDef_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getSortedArrDefComDef_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getSortedArrDefComDef_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getSortedArrDefComDef_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure getSortedArrDefComDef_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getSortedArrDefComDef_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getSortedArrDefComDef_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getSortedArrDefComDef_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getSortedArrDefComDef_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getSortedArrDefComDef_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getSortedArrDefComDef_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getSortedArrDefComDef_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getSortedArrDefComDef_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getSortedArrDefComDef_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getSortedArrDefComDef_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getSortedArrDefComDef_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getSortedArrDefComDef_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getSortedArrDefComDef_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getSortedArrDefComDef_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getSortedArrDefComDef_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getSortedArrDefComDef_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getSortedArrDefComDef_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getSortedArrDefComDef_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getSortedArrDefComDef_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getSortedArrDefComDef_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getSortedArrDefComDef_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getSortedArrDefComDef_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getSortedArrDefComDef_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getSortedArrDefComDef_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getSortedArrDefComDef_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getSortedArrDefComDef_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getSortedArrDefComDef_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getSortedArrDefComDef_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getSortedArrDefComDef_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getSortedArrDefComDef_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ind_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedIndCusComDef_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedIndCusComDef_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedIndCusComDef_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedIndCusComDef_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedIndCusComDef_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedIndCusComDef_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedIndCusComDef_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedIndCusComDef_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedIndCusComDef_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedIndCusComDef_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedIndCusComDef_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedIndCusComDef_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedIndCusComDef_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedIndCusComDef_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedIndCusComDef_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedIndCusComDef_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedIndCusComDef_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedIndCusComDef_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedIndCusComDef_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedIndCusComDef_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedIndCusComDef_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedIndCusComDef_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedIndCusComDef_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedIndCusComDef_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedIndCusComDef_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedIndCusComDef_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedIndCusComDef_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedIndCusComDef_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedIndCusComDef_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedIndCusComDef_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedIndCusComDef_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedIndCusComDef_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedIndCusComDef_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedIndCusComDef_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedIndCusComDef_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedIndCusComDef_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ind_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComDef_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComDef_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComDef_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComDef_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComDef_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrCusComDef_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComDef_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComDef_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComDef_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComDef_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrCusComDef_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrCusComDef_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrCusComDef_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrCusComDef_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrCusComDef_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrCusComDef_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrCusComDef_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrCusComDef_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrCusComDef_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrCusComDef_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrCusComDef_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrCusComDef_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrCusComDef_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrCusComDef_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrCusComDef_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrCusComDef_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrCusComDef_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrCusComDef_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrCusComDef_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrCusComDef_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComDef_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComDef_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComDef_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComDef_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComDef_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrCusComDef_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Qsorti_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComQsorti_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComQsorti_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComQsorti_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComQsorti_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComQsorti_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrCusComQsorti_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComQsorti_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComQsorti_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComQsorti_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComQsorti_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrCusComQsorti_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrCusComQsorti_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrCusComQsorti_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrCusComQsorti_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrCusComQsorti_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrCusComQsorti_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrCusComQsorti_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrCusComQsorti_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrCusComQsorti_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrCusComQsorti_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrCusComQsorti_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrCusComQsorti_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrCusComQsorti_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrCusComQsorti_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrCusComQsorti_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrCusComQsorti_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrCusComQsorti_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrCusComQsorti_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrCusComQsorti_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrCusComQsorti_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComQsorti_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComQsorti_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComQsorti_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComQsorti_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComQsorti_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrCusComQsorti_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Qsorti_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Qsortr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComQsortr_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComQsortr_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComQsortr_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComQsortr_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComQsortr_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrCusComQsortr_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComQsortr_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComQsortr_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComQsortr_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComQsortr_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrCusComQsortr_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrCusComQsortr_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrCusComQsortr_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrCusComQsortr_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrCusComQsortr_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrCusComQsortr_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrCusComQsortr_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrCusComQsortr_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrCusComQsortr_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrCusComQsortr_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrCusComQsortr_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrCusComQsortr_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrCusComQsortr_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrCusComQsortr_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrCusComQsortr_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrCusComQsortr_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrCusComQsortr_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrCusComQsortr_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrCusComQsortr_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrCusComQsortr_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComQsortr_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComQsortr_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComQsortr_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComQsortr_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComQsortr_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrCusComQsortr_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Qsortr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Qsortrdp_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComQsortrdp_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComQsortrdp_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComQsortrdp_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComQsortrdp_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComQsortrdp_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrCusComQsortrdp_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComQsortrdp_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrCusComQsortrdp_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Qsortrdp_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Bubble_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComBubble_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComBubble_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComBubble_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComBubble_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComBubble_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrCusComBubble_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComBubble_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComBubble_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComBubble_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComBubble_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrCusComBubble_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrCusComBubble_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrCusComBubble_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrCusComBubble_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrCusComBubble_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrCusComBubble_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrCusComBubble_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrCusComBubble_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrCusComBubble_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrCusComBubble_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrCusComBubble_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrCusComBubble_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrCusComBubble_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrCusComBubble_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrCusComBubble_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrCusComBubble_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrCusComBubble_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrCusComBubble_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrCusComBubble_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrCusComBubble_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComBubble_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComBubble_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComBubble_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComBubble_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComBubble_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrCusComBubble_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Bubble_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Heapi_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComHeapi_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComHeapi_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComHeapi_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComHeapi_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComHeapi_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrCusComHeapi_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComHeapi_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComHeapi_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComHeapi_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComHeapi_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrCusComHeapi_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrCusComHeapi_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrCusComHeapi_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrCusComHeapi_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrCusComHeapi_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrCusComHeapi_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrCusComHeapi_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrCusComHeapi_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrCusComHeapi_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrCusComHeapi_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrCusComHeapi_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrCusComHeapi_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrCusComHeapi_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrCusComHeapi_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrCusComHeapi_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrCusComHeapi_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrCusComHeapi_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrCusComHeapi_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrCusComHeapi_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrCusComHeapi_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComHeapi_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComHeapi_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComHeapi_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComHeapi_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComHeapi_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrCusComHeapi_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Heapi_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Heapr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComHeapr_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComHeapr_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComHeapr_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComHeapr_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComHeapr_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrCusComHeapr_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComHeapr_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComHeapr_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComHeapr_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComHeapr_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrCusComHeapr_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrCusComHeapr_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrCusComHeapr_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrCusComHeapr_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrCusComHeapr_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrCusComHeapr_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrCusComHeapr_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrCusComHeapr_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrCusComHeapr_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrCusComHeapr_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrCusComHeapr_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrCusComHeapr_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrCusComHeapr_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrCusComHeapr_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrCusComHeapr_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrCusComHeapr_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrCusComHeapr_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrCusComHeapr_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrCusComHeapr_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrCusComHeapr_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComHeapr_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComHeapr_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComHeapr_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComHeapr_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComHeapr_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrCusComHeapr_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Heapr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Insertionl_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComInsertionl_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComInsertionl_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComInsertionl_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComInsertionl_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComInsertionl_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrCusComInsertionl_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComInsertionl_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrCusComInsertionl_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Insertionl_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Insertionb_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComInsertionb_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComInsertionb_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComInsertionb_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComInsertionb_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComInsertionb_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrCusComInsertionb_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComInsertionb_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrCusComInsertionb_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Insertionb_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Merger_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComMerger_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComMerger_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComMerger_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComMerger_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComMerger_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrCusComMerger_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComMerger_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComMerger_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComMerger_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComMerger_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrCusComMerger_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrCusComMerger_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrCusComMerger_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrCusComMerger_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrCusComMerger_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrCusComMerger_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrCusComMerger_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrCusComMerger_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrCusComMerger_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrCusComMerger_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrCusComMerger_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrCusComMerger_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrCusComMerger_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrCusComMerger_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrCusComMerger_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrCusComMerger_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrCusComMerger_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrCusComMerger_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrCusComMerger_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrCusComMerger_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComMerger_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComMerger_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComMerger_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComMerger_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComMerger_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrCusComMerger_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Merger_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Selection_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComSelection_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComSelection_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComSelection_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComSelection_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComSelection_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrCusComSelection_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComSelection_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComSelection_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComSelection_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComSelection_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrCusComSelection_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrCusComSelection_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrCusComSelection_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrCusComSelection_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrCusComSelection_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrCusComSelection_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrCusComSelection_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrCusComSelection_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrCusComSelection_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrCusComSelection_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrCusComSelection_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrCusComSelection_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrCusComSelection_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrCusComSelection_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrCusComSelection_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrCusComSelection_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrCusComSelection_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrCusComSelection_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrCusComSelection_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrCusComSelection_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComSelection_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComSelection_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComSelection_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComSelection_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComSelection_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrCusComSelection_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Selection_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Shell_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComShell_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComShell_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComShell_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComShell_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComShell_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrCusComShell_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComShell_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComShell_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComShell_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComShell_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrCusComShell_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrCusComShell_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrCusComShell_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrCusComShell_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrCusComShell_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrCusComShell_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrCusComShell_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrCusComShell_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrCusComShell_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrCusComShell_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrCusComShell_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrCusComShell_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrCusComShell_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrCusComShell_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrCusComShell_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrCusComShell_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrCusComShell_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrCusComShell_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrCusComShell_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrCusComShell_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrCusComShell_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrCusComShell_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrCusComShell_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrCusComShell_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrCusComShell_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrCusComShell_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Shell_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ind_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedIndDefComDef_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedIndDefComDef_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedIndDefComDef_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedIndDefComDef_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedIndDefComDef_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedIndDefComDef_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedIndDefComDef_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedIndDefComDef_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedIndDefComDef_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedIndDefComDef_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedIndDefComDef_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedIndDefComDef_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedIndDefComDef_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedIndDefComDef_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedIndDefComDef_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedIndDefComDef_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedIndDefComDef_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedIndDefComDef_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedIndDefComDef_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedIndDefComDef_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedIndDefComDef_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedIndDefComDef_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedIndDefComDef_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedIndDefComDef_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedIndDefComDef_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedIndDefComDef_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedIndDefComDef_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedIndDefComDef_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedIndDefComDef_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedIndDefComDef_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedIndDefComDef_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedIndDefComDef_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedIndDefComDef_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedIndDefComDef_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedIndDefComDef_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedIndDefComDef_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ind_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComDef_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComDef_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComDef_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComDef_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComDef_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrDefComDef_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComDef_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComDef_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComDef_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComDef_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrDefComDef_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrDefComDef_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrDefComDef_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrDefComDef_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrDefComDef_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrDefComDef_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrDefComDef_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrDefComDef_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrDefComDef_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrDefComDef_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrDefComDef_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrDefComDef_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrDefComDef_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrDefComDef_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrDefComDef_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrDefComDef_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrDefComDef_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrDefComDef_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrDefComDef_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrDefComDef_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComDef_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComDef_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComDef_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComDef_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComDef_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrDefComDef_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Qsorti_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComQsorti_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComQsorti_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComQsorti_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComQsorti_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComQsorti_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrDefComQsorti_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComQsorti_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComQsorti_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComQsorti_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComQsorti_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrDefComQsorti_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrDefComQsorti_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrDefComQsorti_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrDefComQsorti_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrDefComQsorti_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrDefComQsorti_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrDefComQsorti_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrDefComQsorti_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrDefComQsorti_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrDefComQsorti_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrDefComQsorti_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrDefComQsorti_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrDefComQsorti_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrDefComQsorti_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrDefComQsorti_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrDefComQsorti_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrDefComQsorti_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrDefComQsorti_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrDefComQsorti_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrDefComQsorti_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComQsorti_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComQsorti_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComQsorti_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComQsorti_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComQsorti_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrDefComQsorti_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Qsorti_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Qsortr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComQsortr_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComQsortr_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComQsortr_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComQsortr_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComQsortr_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrDefComQsortr_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComQsortr_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComQsortr_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComQsortr_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComQsortr_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrDefComQsortr_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrDefComQsortr_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrDefComQsortr_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrDefComQsortr_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrDefComQsortr_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrDefComQsortr_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrDefComQsortr_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrDefComQsortr_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrDefComQsortr_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrDefComQsortr_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrDefComQsortr_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrDefComQsortr_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrDefComQsortr_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrDefComQsortr_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrDefComQsortr_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrDefComQsortr_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrDefComQsortr_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrDefComQsortr_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrDefComQsortr_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrDefComQsortr_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComQsortr_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComQsortr_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComQsortr_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComQsortr_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComQsortr_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrDefComQsortr_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Qsortr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Qsortrdp_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComQsortrdp_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComQsortrdp_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComQsortrdp_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComQsortrdp_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComQsortrdp_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrDefComQsortrdp_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComQsortrdp_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrDefComQsortrdp_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Qsortrdp_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Bubble_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComBubble_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComBubble_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComBubble_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComBubble_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComBubble_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrDefComBubble_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComBubble_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComBubble_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComBubble_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComBubble_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrDefComBubble_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrDefComBubble_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrDefComBubble_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrDefComBubble_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrDefComBubble_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrDefComBubble_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrDefComBubble_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrDefComBubble_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrDefComBubble_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrDefComBubble_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrDefComBubble_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrDefComBubble_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrDefComBubble_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrDefComBubble_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrDefComBubble_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrDefComBubble_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrDefComBubble_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrDefComBubble_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrDefComBubble_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrDefComBubble_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComBubble_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComBubble_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComBubble_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComBubble_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComBubble_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrDefComBubble_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Bubble_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Heapi_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComHeapi_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComHeapi_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComHeapi_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComHeapi_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComHeapi_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrDefComHeapi_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComHeapi_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComHeapi_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComHeapi_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComHeapi_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrDefComHeapi_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrDefComHeapi_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrDefComHeapi_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrDefComHeapi_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrDefComHeapi_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrDefComHeapi_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrDefComHeapi_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrDefComHeapi_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrDefComHeapi_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrDefComHeapi_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrDefComHeapi_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrDefComHeapi_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrDefComHeapi_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrDefComHeapi_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrDefComHeapi_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrDefComHeapi_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrDefComHeapi_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrDefComHeapi_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrDefComHeapi_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrDefComHeapi_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComHeapi_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComHeapi_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComHeapi_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComHeapi_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComHeapi_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrDefComHeapi_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Heapi_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Heapr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComHeapr_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComHeapr_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComHeapr_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComHeapr_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComHeapr_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrDefComHeapr_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComHeapr_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComHeapr_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComHeapr_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComHeapr_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrDefComHeapr_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrDefComHeapr_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrDefComHeapr_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrDefComHeapr_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrDefComHeapr_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrDefComHeapr_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrDefComHeapr_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrDefComHeapr_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrDefComHeapr_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrDefComHeapr_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrDefComHeapr_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrDefComHeapr_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrDefComHeapr_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrDefComHeapr_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrDefComHeapr_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrDefComHeapr_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrDefComHeapr_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrDefComHeapr_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrDefComHeapr_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrDefComHeapr_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComHeapr_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComHeapr_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComHeapr_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComHeapr_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComHeapr_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrDefComHeapr_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Heapr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Insertionl_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComInsertionl_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComInsertionl_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComInsertionl_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComInsertionl_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComInsertionl_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrDefComInsertionl_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComInsertionl_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrDefComInsertionl_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Insertionl_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Insertionb_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComInsertionb_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComInsertionb_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComInsertionb_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComInsertionb_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComInsertionb_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrDefComInsertionb_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComInsertionb_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrDefComInsertionb_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Insertionb_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Merger_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComMerger_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComMerger_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComMerger_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComMerger_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComMerger_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrDefComMerger_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComMerger_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComMerger_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComMerger_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComMerger_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrDefComMerger_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrDefComMerger_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrDefComMerger_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrDefComMerger_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrDefComMerger_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrDefComMerger_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrDefComMerger_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrDefComMerger_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrDefComMerger_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrDefComMerger_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrDefComMerger_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrDefComMerger_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrDefComMerger_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrDefComMerger_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrDefComMerger_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrDefComMerger_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrDefComMerger_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrDefComMerger_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrDefComMerger_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrDefComMerger_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComMerger_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComMerger_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComMerger_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComMerger_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComMerger_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrDefComMerger_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Merger_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Selection_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComSelection_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComSelection_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComSelection_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComSelection_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComSelection_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrDefComSelection_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComSelection_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComSelection_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComSelection_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComSelection_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrDefComSelection_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrDefComSelection_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrDefComSelection_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrDefComSelection_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrDefComSelection_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrDefComSelection_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrDefComSelection_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrDefComSelection_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrDefComSelection_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrDefComSelection_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrDefComSelection_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrDefComSelection_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrDefComSelection_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrDefComSelection_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrDefComSelection_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrDefComSelection_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrDefComSelection_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrDefComSelection_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrDefComSelection_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrDefComSelection_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComSelection_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComSelection_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComSelection_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComSelection_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComSelection_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrDefComSelection_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Selection_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Shell_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComShell_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComShell_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComShell_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComShell_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComShell_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
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
    module procedure setSortedArrDefComShell_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComShell_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComShell_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComShell_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComShell_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSortedArrDefComShell_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSortedArrDefComShell_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSortedArrDefComShell_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSortedArrDefComShell_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSortedArrDefComShell_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSortedArrDefComShell_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSortedArrDefComShell_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSortedArrDefComShell_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSortedArrDefComShell_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSortedArrDefComShell_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSortedArrDefComShell_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSortedArrDefComShell_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSortedArrDefComShell_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSortedArrDefComShell_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSortedArrDefComShell_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSortedArrDefComShell_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSortedArrDefComShell_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSortedArrDefComShell_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSortedArrDefComShell_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSortedArrDefComShell_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setSortedArrDefComShell_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSortedArrDefComShell_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSortedArrDefComShell_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSortedArrDefComShell_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSortedArrDefComShell_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setSortedArrDefComShell_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_arraySort@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Shell_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines