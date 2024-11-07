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
!>  This file contains procedure implementations of [pm_matrixInit](@ref pm_matrixInit).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_matrixInit) routines ! LCOV_EXCL_LINE

    !   \bug Bypass Intel ifort bug for too many `use` statements in submodule.
#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif
    use pm_option, only: getOption

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatInit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2XX0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatInitXXD_D2XX0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatInitXXD_D2XX0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatInitXXD_D2XX0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatInitXXD_D2XX0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatInitXXD_D2XX0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatInitXXD_D2XX0_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatInitXXD_D2XX0_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatInitXXD_D2XX0_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatInitXXD_D2XX0_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatInitXXD_D2XX0_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatInitXXD_D2XX0_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatInitXXD_D2XX0_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatInitXXD_D2XX0_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatInitXXD_D2XX0_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatInitXXD_D2XX0_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInitXXD_D2XX0_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInitXXD_D2XX0_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInitXXD_D2XX0_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInitXXD_D2XX0_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInitXXD_D2XX0_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInitXXD_D2XX0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInitXXD_D2XX0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInitXXD_D2XX0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInitXXD_D2XX0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInitXXD_D2XX0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2XX0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2XX1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatInitXXD_D2XX1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatInitXXD_D2XX1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatInitXXD_D2XX1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatInitXXD_D2XX1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatInitXXD_D2XX1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatInitXXD_D2XX1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatInitXXD_D2XX1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatInitXXD_D2XX1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatInitXXD_D2XX1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatInitXXD_D2XX1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatInitXXD_D2XX1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatInitXXD_D2XX1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatInitXXD_D2XX1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatInitXXD_D2XX1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatInitXXD_D2XX1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInitXXD_D2XX1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInitXXD_D2XX1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInitXXD_D2XX1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInitXXD_D2XX1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInitXXD_D2XX1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInitXXD_D2XX1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInitXXD_D2XX1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInitXXD_D2XX1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInitXXD_D2XX1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInitXXD_D2XX1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2XX1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2X00_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatInitXLD_D2X00_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatInitXLD_D2X00_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatInitXLD_D2X00_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatInitXLD_D2X00_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatInitXLD_D2X00_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatInitXLD_D2X00_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatInitXLD_D2X00_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatInitXLD_D2X00_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatInitXLD_D2X00_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatInitXLD_D2X00_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatInitXLD_D2X00_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatInitXLD_D2X00_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatInitXLD_D2X00_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatInitXLD_D2X00_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatInitXLD_D2X00_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInitXLD_D2X00_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInitXLD_D2X00_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInitXLD_D2X00_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInitXLD_D2X00_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInitXLD_D2X00_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInitXLD_D2X00_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInitXLD_D2X00_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInitXLD_D2X00_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInitXLD_D2X00_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInitXLD_D2X00_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2X00_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2X01_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatInitXLD_D2X01_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatInitXLD_D2X01_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatInitXLD_D2X01_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatInitXLD_D2X01_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatInitXLD_D2X01_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatInitXLD_D2X01_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatInitXLD_D2X01_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatInitXLD_D2X01_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatInitXLD_D2X01_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatInitXLD_D2X01_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatInitXLD_D2X01_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatInitXLD_D2X01_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatInitXLD_D2X01_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatInitXLD_D2X01_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatInitXLD_D2X01_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInitXLD_D2X01_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInitXLD_D2X01_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInitXLD_D2X01_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInitXLD_D2X01_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInitXLD_D2X01_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInitXLD_D2X01_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInitXLD_D2X01_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInitXLD_D2X01_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInitXLD_D2X01_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInitXLD_D2X01_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2X01_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D20X0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatInitUXD_D20X0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatInitUXD_D20X0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatInitUXD_D20X0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatInitUXD_D20X0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatInitUXD_D20X0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatInitUXD_D20X0_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatInitUXD_D20X0_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatInitUXD_D20X0_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatInitUXD_D20X0_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatInitUXD_D20X0_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatInitUXD_D20X0_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatInitUXD_D20X0_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatInitUXD_D20X0_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatInitUXD_D20X0_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatInitUXD_D20X0_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInitUXD_D20X0_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInitUXD_D20X0_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInitUXD_D20X0_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInitUXD_D20X0_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInitUXD_D20X0_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInitUXD_D20X0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInitUXD_D20X0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInitUXD_D20X0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInitUXD_D20X0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInitUXD_D20X0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D20X0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D20X1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatInitUXD_D20X1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatInitUXD_D20X1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatInitUXD_D20X1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatInitUXD_D20X1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatInitUXD_D20X1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatInitUXD_D20X1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatInitUXD_D20X1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatInitUXD_D20X1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatInitUXD_D20X1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatInitUXD_D20X1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatInitUXD_D20X1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatInitUXD_D20X1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatInitUXD_D20X1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatInitUXD_D20X1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatInitUXD_D20X1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInitUXD_D20X1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInitUXD_D20X1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInitUXD_D20X1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInitUXD_D20X1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInitUXD_D20X1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInitUXD_D20X1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInitUXD_D20X1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInitUXD_D20X1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInitUXD_D20X1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInitUXD_D20X1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D20X1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D200X_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatInitULX_D200X_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatInitULX_D200X_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatInitULX_D200X_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatInitULX_D200X_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatInitULX_D200X_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatInitULX_D200X_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatInitULX_D200X_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatInitULX_D200X_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatInitULX_D200X_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatInitULX_D200X_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatInitULX_D200X_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatInitULX_D200X_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatInitULX_D200X_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatInitULX_D200X_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatInitULX_D200X_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInitULX_D200X_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInitULX_D200X_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInitULX_D200X_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInitULX_D200X_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInitULX_D200X_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInitULX_D200X_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInitULX_D200X_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInitULX_D200X_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInitULX_D200X_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInitULX_D200X_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D200X_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2000_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatInitULD_D2000_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatInitULD_D2000_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatInitULD_D2000_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatInitULD_D2000_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatInitULD_D2000_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatInitULD_D2000_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatInitULD_D2000_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatInitULD_D2000_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatInitULD_D2000_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatInitULD_D2000_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatInitULD_D2000_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatInitULD_D2000_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatInitULD_D2000_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatInitULD_D2000_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatInitULD_D2000_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInitULD_D2000_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInitULD_D2000_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInitULD_D2000_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInitULD_D2000_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInitULD_D2000_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInitULD_D2000_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInitULD_D2000_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInitULD_D2000_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInitULD_D2000_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInitULD_D2000_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2000_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2001_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMatInitULD_D2001_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMatInitULD_D2001_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMatInitULD_D2001_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMatInitULD_D2001_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMatInitULD_D2001_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatInitULD_D2001_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatInitULD_D2001_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatInitULD_D2001_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatInitULD_D2001_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatInitULD_D2001_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMatInitULD_D2001_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMatInitULD_D2001_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMatInitULD_D2001_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMatInitULD_D2001_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMatInitULD_D2001_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInitULD_D2001_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInitULD_D2001_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInitULD_D2001_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInitULD_D2001_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInitULD_D2001_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInitULD_D2001_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInitULD_D2001_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInitULD_D2001_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInitULD_D2001_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInitULD_D2001_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2001_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMatInit_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatInit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define EXP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2X0X_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_EXP_XLX_D2X0X_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2X0X_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLX_ENABLED

#define UXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D20XX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_EXP_UXX_D20XX_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D20XX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2XXF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_EXP_XXD_D2XXF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2XXF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2XX0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2XX0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2XX1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_EXP_XXD_D2XX1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2XX1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2X00_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_EXP_XLD_D2X00_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2X00_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2X01_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_EXP_XLD_D2X01_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2X01_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D20X0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_EXP_UXD_D20X0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D20X0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D20X1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_EXP_UXD_D20X1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D20X1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D200X_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_EXP_ULX_D200X_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D200X_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2000_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_EXP_ULD_D2000_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2000_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2001_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_EXP_ULD_D2001_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2001_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef EXP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IMP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2X0X_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_IMP_XLX_D2X0X_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2X0X_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLX_ENABLED

#define UXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D20XX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_IMP_UXX_D20XX_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D20XX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2XXF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_IMP_XXD_D2XXF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2XXF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2XX0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2XX0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2XX1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_IMP_XXD_D2XX1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2XX1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2X00_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_IMP_XLD_D2X00_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2X00_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2X01_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_IMP_XLD_D2X01_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2X01_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D20X0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_IMP_UXD_D20X0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D20X0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D20X1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_IMP_UXD_D20X1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D20X1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D200X_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_IMP_ULX_D200X_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D200X_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2000_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_IMP_ULD_D2000_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2000_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2001_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_SK5
        use pm_kind, only: SKG => SK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_SK4
        use pm_kind, only: SKG => SK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_SK3
        use pm_kind, only: SKG => SK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_SK2
        use pm_kind, only: SKG => SK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_SK1
        use pm_kind, only: SKG => SK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_LK5
        use pm_kind, only: LKG => LK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_LK4
        use pm_kind, only: LKG => LK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_LK3
        use pm_kind, only: LKG => LK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_LK2
        use pm_kind, only: LKG => LK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_LK1
        use pm_kind, only: LKG => LK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInit_IMP_ULD_D2001_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2001_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IMP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatInit_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  CHECK_ASSERTION

end submodule routines
