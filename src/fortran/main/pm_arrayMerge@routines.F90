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
!>  This file contains procedure implementations of [pm_arrayMerge](@ref pm_arrayMerge).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_arrayMerge) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
    use pm_arraySort, only: isAscending
    use pm_arraySort, only: isSortedCheck => isSorted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMerged_ENABLED 1

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
    module procedure getMergedDefCom_D0_SK5
        use pm_kind, only: TKG => SK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMergedDefCom_D0_SK4
        use pm_kind, only: TKG => SK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMergedDefCom_D0_SK3
        use pm_kind, only: TKG => SK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMergedDefCom_D0_SK2
        use pm_kind, only: TKG => SK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMergedDefCom_D0_SK1
        use pm_kind, only: TKG => SK1
#include "pm_arrayMerge@routines.inc.F90"
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
    module procedure getMergedDefCom_D1_SK5
        use pm_kind, only: TKG => SK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMergedDefCom_D1_SK4
        use pm_kind, only: TKG => SK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMergedDefCom_D1_SK3
        use pm_kind, only: TKG => SK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMergedDefCom_D1_SK2
        use pm_kind, only: TKG => SK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMergedDefCom_D1_SK1
        use pm_kind, only: TKG => SK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMergedDefCom_D1_IK5
        use pm_kind, only: TKG => IK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMergedDefCom_D1_IK4
        use pm_kind, only: TKG => IK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMergedDefCom_D1_IK3
        use pm_kind, only: TKG => IK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMergedDefCom_D1_IK2
        use pm_kind, only: TKG => IK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMergedDefCom_D1_IK1
        use pm_kind, only: TKG => IK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMergedDefCom_D1_LK5
        use pm_kind, only: TKG => LK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMergedDefCom_D1_LK4
        use pm_kind, only: TKG => LK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMergedDefCom_D1_LK3
        use pm_kind, only: TKG => LK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMergedDefCom_D1_LK2
        use pm_kind, only: TKG => LK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMergedDefCom_D1_LK1
        use pm_kind, only: TKG => LK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMergedDefCom_D1_CK5
        use pm_kind, only: TKG => CK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMergedDefCom_D1_CK4
        use pm_kind, only: TKG => CK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMergedDefCom_D1_CK3
        use pm_kind, only: TKG => CK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMergedDefCom_D1_CK2
        use pm_kind, only: TKG => CK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMergedDefCom_D1_CK1
        use pm_kind, only: TKG => CK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMergedDefCom_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMergedDefCom_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMergedDefCom_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMergedDefCom_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMergedDefCom_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! LCOV_EXCL_START
#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getMergedDefCom_D1_PSSK5
        use pm_kind, only: TKG => SK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMergedDefCom_D1_PSSK4
        use pm_kind, only: TKG => SK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMergedDefCom_D1_PSSK3
        use pm_kind, only: TKG => SK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMergedDefCom_D1_PSSK2
        use pm_kind, only: TKG => SK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMergedDefCom_D1_PSSK1
        use pm_kind, only: TKG => SK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

! LCOV_EXCL_STOP

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getMergedDefCom_D1_BSSK
        use pm_kind, only: TKG => SK
#include "pm_arrayMerge@routines.inc.F90"
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
    module procedure getMergedCusCom_D0_SK5
        use pm_kind, only: TKG => SK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMergedCusCom_D0_SK4
        use pm_kind, only: TKG => SK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMergedCusCom_D0_SK3
        use pm_kind, only: TKG => SK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMergedCusCom_D0_SK2
        use pm_kind, only: TKG => SK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMergedCusCom_D0_SK1
        use pm_kind, only: TKG => SK1
#include "pm_arrayMerge@routines.inc.F90"
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
    module procedure getMergedCusCom_D1_SK5
        use pm_kind, only: TKG => SK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMergedCusCom_D1_SK4
        use pm_kind, only: TKG => SK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMergedCusCom_D1_SK3
        use pm_kind, only: TKG => SK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMergedCusCom_D1_SK2
        use pm_kind, only: TKG => SK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMergedCusCom_D1_SK1
        use pm_kind, only: TKG => SK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMergedCusCom_D1_IK5
        use pm_kind, only: TKG => IK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMergedCusCom_D1_IK4
        use pm_kind, only: TKG => IK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMergedCusCom_D1_IK3
        use pm_kind, only: TKG => IK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMergedCusCom_D1_IK2
        use pm_kind, only: TKG => IK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMergedCusCom_D1_IK1
        use pm_kind, only: TKG => IK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMergedCusCom_D1_LK5
        use pm_kind, only: TKG => LK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMergedCusCom_D1_LK4
        use pm_kind, only: TKG => LK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMergedCusCom_D1_LK3
        use pm_kind, only: TKG => LK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMergedCusCom_D1_LK2
        use pm_kind, only: TKG => LK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMergedCusCom_D1_LK1
        use pm_kind, only: TKG => LK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMergedCusCom_D1_CK5
        use pm_kind, only: TKG => CK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMergedCusCom_D1_CK4
        use pm_kind, only: TKG => CK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMergedCusCom_D1_CK3
        use pm_kind, only: TKG => CK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMergedCusCom_D1_CK2
        use pm_kind, only: TKG => CK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMergedCusCom_D1_CK1
        use pm_kind, only: TKG => CK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMergedCusCom_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMergedCusCom_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMergedCusCom_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMergedCusCom_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMergedCusCom_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! LCOV_EXCL_START

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getMergedCusCom_D1_PSSK5
        use pm_kind, only: TKG => SK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMergedCusCom_D1_PSSK4
        use pm_kind, only: TKG => SK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMergedCusCom_D1_PSSK3
        use pm_kind, only: TKG => SK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMergedCusCom_D1_PSSK2
        use pm_kind, only: TKG => SK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMergedCusCom_D1_PSSK1
        use pm_kind, only: TKG => SK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

! LCOV_EXCL_STOP

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getMergedCusCom_D1_BSSK
        use pm_kind, only: TKG => SK
#include "pm_arrayMerge@routines.inc.F90"
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

#undef getMerged_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMerged_ENABLED 1

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
    module procedure setMergedDefCom_D0_SK5
        use pm_kind, only: TKG => SK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMergedDefCom_D0_SK4
        use pm_kind, only: TKG => SK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMergedDefCom_D0_SK3
        use pm_kind, only: TKG => SK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMergedDefCom_D0_SK2
        use pm_kind, only: TKG => SK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMergedDefCom_D0_SK1
        use pm_kind, only: TKG => SK1
#include "pm_arrayMerge@routines.inc.F90"
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
    module procedure setMergedDefCom_D1_SK5
        use pm_kind, only: TKG => SK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMergedDefCom_D1_SK4
        use pm_kind, only: TKG => SK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMergedDefCom_D1_SK3
        use pm_kind, only: TKG => SK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMergedDefCom_D1_SK2
        use pm_kind, only: TKG => SK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMergedDefCom_D1_SK1
        use pm_kind, only: TKG => SK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMergedDefCom_D1_IK5
        use pm_kind, only: TKG => IK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMergedDefCom_D1_IK4
        use pm_kind, only: TKG => IK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMergedDefCom_D1_IK3
        use pm_kind, only: TKG => IK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMergedDefCom_D1_IK2
        use pm_kind, only: TKG => IK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMergedDefCom_D1_IK1
        use pm_kind, only: TKG => IK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMergedDefCom_D1_LK5
        use pm_kind, only: TKG => LK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMergedDefCom_D1_LK4
        use pm_kind, only: TKG => LK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMergedDefCom_D1_LK3
        use pm_kind, only: TKG => LK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMergedDefCom_D1_LK2
        use pm_kind, only: TKG => LK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMergedDefCom_D1_LK1
        use pm_kind, only: TKG => LK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMergedDefCom_D1_CK5
        use pm_kind, only: TKG => CK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMergedDefCom_D1_CK4
        use pm_kind, only: TKG => CK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMergedDefCom_D1_CK3
        use pm_kind, only: TKG => CK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMergedDefCom_D1_CK2
        use pm_kind, only: TKG => CK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMergedDefCom_D1_CK1
        use pm_kind, only: TKG => CK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMergedDefCom_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMergedDefCom_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMergedDefCom_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMergedDefCom_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMergedDefCom_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! LCOV_EXCL_START

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setMergedDefCom_D1_PSSK5
        use pm_kind, only: TKG => SK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMergedDefCom_D1_PSSK4
        use pm_kind, only: TKG => SK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMergedDefCom_D1_PSSK3
        use pm_kind, only: TKG => SK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMergedDefCom_D1_PSSK2
        use pm_kind, only: TKG => SK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMergedDefCom_D1_PSSK1
        use pm_kind, only: TKG => SK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

! LCOV_EXCL_STOP

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setMergedDefCom_D1_BSSK
        use pm_kind, only: TKG => SK
#include "pm_arrayMerge@routines.inc.F90"
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
    module procedure setMergedCusCom_D0_SK5
        use pm_kind, only: TKG => SK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMergedCusCom_D0_SK4
        use pm_kind, only: TKG => SK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMergedCusCom_D0_SK3
        use pm_kind, only: TKG => SK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMergedCusCom_D0_SK2
        use pm_kind, only: TKG => SK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMergedCusCom_D0_SK1
        use pm_kind, only: TKG => SK1
#include "pm_arrayMerge@routines.inc.F90"
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
    module procedure setMergedCusCom_D1_SK5
        use pm_kind, only: TKG => SK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMergedCusCom_D1_SK4
        use pm_kind, only: TKG => SK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMergedCusCom_D1_SK3
        use pm_kind, only: TKG => SK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMergedCusCom_D1_SK2
        use pm_kind, only: TKG => SK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMergedCusCom_D1_SK1
        use pm_kind, only: TKG => SK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMergedCusCom_D1_IK5
        use pm_kind, only: TKG => IK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMergedCusCom_D1_IK4
        use pm_kind, only: TKG => IK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMergedCusCom_D1_IK3
        use pm_kind, only: TKG => IK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMergedCusCom_D1_IK2
        use pm_kind, only: TKG => IK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMergedCusCom_D1_IK1
        use pm_kind, only: TKG => IK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMergedCusCom_D1_LK5
        use pm_kind, only: TKG => LK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMergedCusCom_D1_LK4
        use pm_kind, only: TKG => LK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMergedCusCom_D1_LK3
        use pm_kind, only: TKG => LK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMergedCusCom_D1_LK2
        use pm_kind, only: TKG => LK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMergedCusCom_D1_LK1
        use pm_kind, only: TKG => LK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMergedCusCom_D1_CK5
        use pm_kind, only: TKG => CK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMergedCusCom_D1_CK4
        use pm_kind, only: TKG => CK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMergedCusCom_D1_CK3
        use pm_kind, only: TKG => CK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMergedCusCom_D1_CK2
        use pm_kind, only: TKG => CK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMergedCusCom_D1_CK1
        use pm_kind, only: TKG => CK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMergedCusCom_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMergedCusCom_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMergedCusCom_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMergedCusCom_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMergedCusCom_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! LCOV_EXCL_START

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setMergedCusCom_D1_PSSK5
        use pm_kind, only: TKG => SK5
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMergedCusCom_D1_PSSK4
        use pm_kind, only: TKG => SK4
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMergedCusCom_D1_PSSK3
        use pm_kind, only: TKG => SK3
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMergedCusCom_D1_PSSK2
        use pm_kind, only: TKG => SK2
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMergedCusCom_D1_PSSK1
        use pm_kind, only: TKG => SK1
#include "pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

! LCOV_EXCL_STOP

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setMergedCusCom_D1_BSSK
        use pm_kind, only: TKG => SK
#include "pm_arrayMerge@routines.inc.F90"
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

#undef setMerged_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines
