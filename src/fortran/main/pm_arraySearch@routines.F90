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
!>  This file contains procedure implementations of [pm_arraySearch](@ref pm_arraySearch).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_arraySearch) routines ! LCOV_EXCL_LINE

    use pm_kind, only: LK

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

#define getBin_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getBinDefCom_D0_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getBinDefCom_D0_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getBinDefCom_D0_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getBinDefCom_D0_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getBinDefCom_D0_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getBinDefCom_D1_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getBinDefCom_D1_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getBinDefCom_D1_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getBinDefCom_D1_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getBinDefCom_D1_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getBinDefCom_D1_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getBinDefCom_D1_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getBinDefCom_D1_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getBinDefCom_D1_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getBinDefCom_D1_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getBinDefCom_D1_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getBinDefCom_D1_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getBinDefCom_D1_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getBinDefCom_D1_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getBinDefCom_D1_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getBinDefCom_D1_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getBinDefCom_D1_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getBinDefCom_D1_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getBinDefCom_D1_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getBinDefCom_D1_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CusCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getBinCusCom_D0_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getBinCusCom_D0_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getBinCusCom_D0_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getBinCusCom_D0_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getBinCusCom_D0_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getBinCusCom_D1_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getBinCusCom_D1_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getBinCusCom_D1_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getBinCusCom_D1_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getBinCusCom_D1_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getBinCusCom_D1_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getBinCusCom_D1_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getBinCusCom_D1_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getBinCusCom_D1_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getBinCusCom_D1_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getBinCusCom_D1_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getBinCusCom_D1_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getBinCusCom_D1_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getBinCusCom_D1_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getBinCusCom_D1_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getBinCusCom_D1_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getBinCusCom_D1_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getBinCusCom_D1_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getBinCusCom_D1_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getBinCusCom_D1_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_arraySearch@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getBin_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines