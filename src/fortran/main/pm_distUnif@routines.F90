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
!>  This file contains procedure implementations of [pm_distUnif](@ref pm_distUnif).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 PM, September 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distUnif) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_kind, only: RKB
    implicit none

    !integer(IK) :: xoshiro256ssStreamLenMinusOne = xoshiro256ssStreamBitSize - 1_IK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getUnifCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifCDF_DD_IK5
        use pm_kind, only: RKC => RK, IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifCDF_DD_IK4
        use pm_kind, only: RKC => RK, IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifCDF_DD_IK3
        use pm_kind, only: RKC => RK, IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifCDF_DD_IK2
        use pm_kind, only: RKC => RK, IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifCDF_DD_IK1
        use pm_kind, only: RKC => RK, IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifCDF_DD_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifCDF_DD_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifCDF_DD_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifCDF_DD_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifCDF_DD_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifCDF_DD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifCDF_DD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifCDF_DD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifCDF_DD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifCDF_DD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifCDF_LU_IK5
        use pm_kind, only: RKC => RK, IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifCDF_LU_IK4
        use pm_kind, only: RKC => RK, IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifCDF_LU_IK3
        use pm_kind, only: RKC => RK, IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifCDF_LU_IK2
        use pm_kind, only: RKC => RK, IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifCDF_LU_IK1
        use pm_kind, only: RKC => RK, IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifCDF_LU_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifCDF_LU_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifCDF_LU_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifCDF_LU_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifCDF_LU_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifCDF_LU_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifCDF_LU_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifCDF_LU_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifCDF_LU_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifCDF_LU_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getUnifCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setUnifCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#if RK5_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_DD_D0_RK5_IK5
        use pm_kind, only: IKC => IK5, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_DD_D0_RK4_IK5
        use pm_kind, only: IKC => IK5, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_DD_D0_RK3_IK5
        use pm_kind, only: IKC => IK5, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_DD_D0_RK2_IK5
        use pm_kind, only: IKC => IK5, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_DD_D0_RK1_IK5
        use pm_kind, only: IKC => IK5, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_DD_D0_RK5_IK4
        use pm_kind, only: IKC => IK4, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_DD_D0_RK4_IK4
        use pm_kind, only: IKC => IK4, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_DD_D0_RK3_IK4
        use pm_kind, only: IKC => IK4, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_DD_D0_RK2_IK4
        use pm_kind, only: IKC => IK4, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_DD_D0_RK1_IK4
        use pm_kind, only: IKC => IK4, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_DD_D0_RK5_IK3
        use pm_kind, only: IKC => IK3, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_DD_D0_RK4_IK3
        use pm_kind, only: IKC => IK3, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_DD_D0_RK3_IK3
        use pm_kind, only: IKC => IK3, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_DD_D0_RK2_IK3
        use pm_kind, only: IKC => IK3, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_DD_D0_RK1_IK3
        use pm_kind, only: IKC => IK3, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_DD_D0_RK5_IK2
        use pm_kind, only: IKC => IK2, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_DD_D0_RK4_IK2
        use pm_kind, only: IKC => IK2, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_DD_D0_RK3_IK2
        use pm_kind, only: IKC => IK2, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_DD_D0_RK2_IK2
        use pm_kind, only: IKC => IK2, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_DD_D0_RK1_IK2
        use pm_kind, only: IKC => IK2, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_DD_D0_RK5_IK1
        use pm_kind, only: IKC => IK1, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_DD_D0_RK4_IK1
        use pm_kind, only: IKC => IK1, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_DD_D0_RK3_IK1
        use pm_kind, only: IKC => IK1, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_DD_D0_RK2_IK1
        use pm_kind, only: IKC => IK1, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_DD_D0_RK1_IK1
        use pm_kind, only: IKC => IK1, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module procedure setUnifCDF_DD_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifCDF_DD_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifCDF_DD_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifCDF_DD_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifCDF_DD_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module procedure setUnifCDF_DD_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifCDF_DD_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifCDF_DD_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifCDF_DD_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifCDF_DD_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_DD_D1_RK5_IK5
        use pm_kind, only: IKC => IK5, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_DD_D1_RK4_IK5
        use pm_kind, only: IKC => IK5, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_DD_D1_RK3_IK5
        use pm_kind, only: IKC => IK5, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_DD_D1_RK2_IK5
        use pm_kind, only: IKC => IK5, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_DD_D1_RK1_IK5
        use pm_kind, only: IKC => IK5, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_DD_D1_RK5_IK4
        use pm_kind, only: IKC => IK4, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_DD_D1_RK4_IK4
        use pm_kind, only: IKC => IK4, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_DD_D1_RK3_IK4
        use pm_kind, only: IKC => IK4, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_DD_D1_RK2_IK4
        use pm_kind, only: IKC => IK4, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_DD_D1_RK1_IK4
        use pm_kind, only: IKC => IK4, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_DD_D1_RK5_IK3
        use pm_kind, only: IKC => IK3, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_DD_D1_RK4_IK3
        use pm_kind, only: IKC => IK3, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_DD_D1_RK3_IK3
        use pm_kind, only: IKC => IK3, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_DD_D1_RK2_IK3
        use pm_kind, only: IKC => IK3, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_DD_D1_RK1_IK3
        use pm_kind, only: IKC => IK3, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_DD_D1_RK5_IK2
        use pm_kind, only: IKC => IK2, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_DD_D1_RK4_IK2
        use pm_kind, only: IKC => IK2, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_DD_D1_RK3_IK2
        use pm_kind, only: IKC => IK2, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_DD_D1_RK2_IK2
        use pm_kind, only: IKC => IK2, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_DD_D1_RK1_IK2
        use pm_kind, only: IKC => IK2, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_DD_D1_RK5_IK1
        use pm_kind, only: IKC => IK1, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_DD_D1_RK4_IK1
        use pm_kind, only: IKC => IK1, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_DD_D1_RK3_IK1
        use pm_kind, only: IKC => IK1, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_DD_D1_RK2_IK1
        use pm_kind, only: IKC => IK1, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_DD_D1_RK1_IK1
        use pm_kind, only: IKC => IK1, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module procedure setUnifCDF_DD_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifCDF_DD_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifCDF_DD_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifCDF_DD_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifCDF_DD_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module procedure setUnifCDF_DD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifCDF_DD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifCDF_DD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifCDF_DD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifCDF_DD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_LU_D0_RK5_IK5
        use pm_kind, only: IKC => IK5, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_LU_D0_RK4_IK5
        use pm_kind, only: IKC => IK5, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_LU_D0_RK3_IK5
        use pm_kind, only: IKC => IK5, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_LU_D0_RK2_IK5
        use pm_kind, only: IKC => IK5, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_LU_D0_RK1_IK5
        use pm_kind, only: IKC => IK5, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_LU_D0_RK5_IK4
        use pm_kind, only: IKC => IK4, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_LU_D0_RK4_IK4
        use pm_kind, only: IKC => IK4, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_LU_D0_RK3_IK4
        use pm_kind, only: IKC => IK4, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_LU_D0_RK2_IK4
        use pm_kind, only: IKC => IK4, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_LU_D0_RK1_IK4
        use pm_kind, only: IKC => IK4, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_LU_D0_RK5_IK3
        use pm_kind, only: IKC => IK3, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_LU_D0_RK4_IK3
        use pm_kind, only: IKC => IK3, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_LU_D0_RK3_IK3
        use pm_kind, only: IKC => IK3, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_LU_D0_RK2_IK3
        use pm_kind, only: IKC => IK3, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_LU_D0_RK1_IK3
        use pm_kind, only: IKC => IK3, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_LU_D0_RK5_IK2
        use pm_kind, only: IKC => IK2, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_LU_D0_RK4_IK2
        use pm_kind, only: IKC => IK2, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_LU_D0_RK3_IK2
        use pm_kind, only: IKC => IK2, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_LU_D0_RK2_IK2
        use pm_kind, only: IKC => IK2, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_LU_D0_RK1_IK2
        use pm_kind, only: IKC => IK2, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_LU_D0_RK5_IK1
        use pm_kind, only: IKC => IK1, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_LU_D0_RK4_IK1
        use pm_kind, only: IKC => IK1, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_LU_D0_RK3_IK1
        use pm_kind, only: IKC => IK1, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_LU_D0_RK2_IK1
        use pm_kind, only: IKC => IK1, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_LU_D0_RK1_IK1
        use pm_kind, only: IKC => IK1, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module procedure setUnifCDF_LU_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifCDF_LU_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifCDF_LU_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifCDF_LU_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifCDF_LU_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module procedure setUnifCDF_LU_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifCDF_LU_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifCDF_LU_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifCDF_LU_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifCDF_LU_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_LU_D1_RK5_IK5
        use pm_kind, only: IKC => IK5, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_LU_D1_RK4_IK5
        use pm_kind, only: IKC => IK5, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_LU_D1_RK3_IK5
        use pm_kind, only: IKC => IK5, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_LU_D1_RK2_IK5
        use pm_kind, only: IKC => IK5, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK5_ENABLED
    module procedure setUnifCDF_LU_D1_RK1_IK5
        use pm_kind, only: IKC => IK5, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_LU_D1_RK5_IK4
        use pm_kind, only: IKC => IK4, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_LU_D1_RK4_IK4
        use pm_kind, only: IKC => IK4, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_LU_D1_RK3_IK4
        use pm_kind, only: IKC => IK4, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_LU_D1_RK2_IK4
        use pm_kind, only: IKC => IK4, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK4_ENABLED
    module procedure setUnifCDF_LU_D1_RK1_IK4
        use pm_kind, only: IKC => IK4, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_LU_D1_RK5_IK3
        use pm_kind, only: IKC => IK3, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_LU_D1_RK4_IK3
        use pm_kind, only: IKC => IK3, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_LU_D1_RK3_IK3
        use pm_kind, only: IKC => IK3, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_LU_D1_RK2_IK3
        use pm_kind, only: IKC => IK3, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK3_ENABLED
    module procedure setUnifCDF_LU_D1_RK1_IK3
        use pm_kind, only: IKC => IK3, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_LU_D1_RK5_IK2
        use pm_kind, only: IKC => IK2, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_LU_D1_RK4_IK2
        use pm_kind, only: IKC => IK2, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_LU_D1_RK3_IK2
        use pm_kind, only: IKC => IK2, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_LU_D1_RK2_IK2
        use pm_kind, only: IKC => IK2, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK2_ENABLED
    module procedure setUnifCDF_LU_D1_RK1_IK2
        use pm_kind, only: IKC => IK2, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_LU_D1_RK5_IK1
        use pm_kind, only: IKC => IK1, RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_LU_D1_RK4_IK1
        use pm_kind, only: IKC => IK1, RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_LU_D1_RK3_IK1
        use pm_kind, only: IKC => IK1, RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_LU_D1_RK2_IK1
        use pm_kind, only: IKC => IK1, RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK1_ENABLED
    module procedure setUnifCDF_LU_D1_RK1_IK1
        use pm_kind, only: IKC => IK1, RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module procedure setUnifCDF_LU_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifCDF_LU_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifCDF_LU_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifCDF_LU_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifCDF_LU_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module procedure setUnifCDF_LU_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifCDF_LU_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifCDF_LU_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifCDF_LU_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifCDF_LU_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setUnifCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define splitmix64_typer_ENABLED 1
    module procedure splitmix64_typer
#include "pm_distUnif@routines.inc.F90"
    end procedure
#undef splitmix64_typer_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define xoshiro256ssg_typer_ENABLED 1
    module procedure xoshiro256ssg_typer
#include "pm_distUnif@routines.inc.F90"
    end procedure
#undef xoshiro256ssg_typer_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define xoshiro256ssw_typer_ENABLED 1
    module procedure xoshiro256ssw_typer
#include "pm_distUnif@routines.inc.F90"
    end procedure
#undef xoshiro256ssw_typer_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setStateNext_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SM64_ENABLED 1
    module procedure setStateNextSM64
#include "pm_distUnif@routines.inc.F90"
    end procedure
#undef SM64_ENABLED

#define X256SSG_ENABLED 1
    module procedure setStateNextX256SSG
#include "pm_distUnif@routines.inc.F90"
    end procedure
#undef X256SSG_ENABLED

#define X256SSW_ENABLED 1
    module procedure setStateNextX256SSW
#include "pm_distUnif@routines.inc.F90"
    end procedure
#undef X256SSW_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setStateNext_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setStateJump_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define X256SSG_ENABLED 1

#define DJ_ENABLED 1
    module procedure setStateJumpX256SSGDJ
#include "pm_distUnif@routines.inc.F90"
    end procedure
#undef DJ_ENABLED

#define AJ_ENABLED 1
    module procedure setStateJumpX256SSGAJ
#include "pm_distUnif@routines.inc.F90"
    end procedure
#undef AJ_ENABLED

#undef X256SSG_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define X256SSW_ENABLED 1

#define DJ_ENABLED 1
    module procedure setStateJumpX256SSWDJ
#include "pm_distUnif@routines.inc.F90"
    end procedure
#undef DJ_ENABLED

#define AJ_ENABLED 1
    module procedure setStateJumpX256SSWAJ
#include "pm_distUnif@routines.inc.F90"
    end procedure
#undef AJ_ENABLED

#undef X256SSW_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setStateJump_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure rngf_typer
        call setUnifRandState(seed, imageID)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getUnifRandStateSizeDef
        call random_seed(size = unifRandStateSize)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getUnifRandStateDef
        allocate(unifRandState(getUnifRandStateSize()))
        if (present(seed) .or. present(imageID)) call setUnifRandState(seed, imageID)
        call random_seed(get = unifRandState)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure setUnifRandStateDef

        if (present(seed) .or. present(imageID)) then

            block

                integer, parameter :: LARGE = huge(0)
                integer, allocatable :: unifRandState(:)
                integer :: offsetImageunifRandState, iseed
                integer :: def_seed, unifRandStateSize
                integer(IK) :: seed_def

                call random_seed(size = unifRandStateSize)
                allocate(unifRandState(0 : unifRandStateSize - 1))

                if (present(seed)) then
                    seed_def = abs(seed)
                else ! Initialize the seed to a random value.
                    block
                        integer(IK) :: values(8)
                        call date_and_time(values = values)
                        seed_def = abs(sum(values))
                    end block
                end if

                ! ensure no overflow occurs in passing to intrinsic functions.
                if (huge(0) < huge(0_IK)) then
                    do
                        if (seed_def < huge(0)) exit
                        seed_def = seed_def - huge(0)
                    end do
                end if

                if (present(imageID)) then
                offsetImageunifRandState = 127 * unifRandStateSize * (int(imageID, IK) - 1)
                else
                    offsetImageunifRandState = 0
                end if

                ! Use seed to construct the random seed on all images

                def_seed = int(seed_def, IK)
                do iseed = 0, unifRandStateSize - 1
                    unifRandState(iseed) = LARGE - def_seed - offsetImageunifRandState - 127 * iseed
                    if (unifRandState(iseed) < 0) then
                        unifRandState(iseed) = -unifRandState(iseed)
                    else
                        unifRandState(iseed) = LARGE - unifRandState(iseed)
                    end if
                end do
                call random_seed(put = unifRandState)
                deallocate(unifRandState)

            end block

            !block
            !call random_init(repeatable = .false._LK, image_distinct = .false._LK)
            !write(*,"(*(g0,:,' '))")
            !write(*,"(*(g0,:,' '))") "unifRandState%val", unifRandState%val
            !write(*,"(*(g0,:,' '))")
            !end block
            !
            ! ATTN: xxx Intel compilers - for some unknown reason, the first generated random number seems to be garbage
            ! so here, the random number generator is iterated a couple of times before further usage.
            ! This needs to be taken care of, in the future. This problem showed itself when StartPoint in ParaDRAM sampler were to be set randomly.
            ! This is where the first instance of random number usage occurs in ParaDRAM sampler.
            ! write(*,*) "unifRandStateObj%imageID, co_unifRandState(1)%val(:): ", unifRandStateObj%imageID, co_unifRandState(1)%val(:)
            !
            ! ATTN: A follow-up on the above issue with the Intel compiler which seems to be a compiler bug: In a truly bizarre behavior,
            ! the Intel compiler random numbers as generated by call random_number() in the pm_statistics module, for example when called from
            ! ParaDRAM_proposal_pmod.inc.F90, are not repeatable even after resetting the random_seed. Even more bizarre is the observation that the
            ! repeatability of the random numbers depends on the loop length (for example as implemented in the debugging of [DistNorm_pmod::getNormRand()](@ref pm_distNorm::getNormRand()).
            ! The same behavior is also observed below, where any loop length less than ~30 yields non-repeatable random number sequences.
            ! This needs an in-depth investigation. Update: Such behavior was also observed with the GNU compiler.
            ! 101 is the number that fixes this issue for both compilers.
            block
                real :: unifrnd(101)
                call random_number(unifrnd)
                !block
                !integer(IK), allocatable :: unifRandStateValue(:)
                !allocate(unifRandStateValue(unifRandState%size))
                !call random_seed(get=unifRandStateValue)
                !write(*,"(*(g0,:,' '))") "unifrnd", unifrnd, unifRandStateValue
                !end block
                !if (this_image()==1) then
                !    write(*,*) "unifRandStateObj%imageID, unifrnd: ", unifrnd
                !    sync images(*)
                !else
                !    sync images(1)
                !    write(*,*) "unifRandStateObj%imageID, unifrnd: ", unifrnd
                !end if
                !if (this_image()==1) read(*,*)
                !sync all
            end block

        else

            call random_seed()

        end if

    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getUnifRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1
    module procedure getUnifRandRNGDDD_D0_LK
        use pm_kind, only: LKC => LK
#include "pm_distUnif@routines.inc.F90"
    end procedure
#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGDLU_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGDLU_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGDLU_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGDLU_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGDLU_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGDLU_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGDLU_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGDLU_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGDLU_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGDLU_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGDLU_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGDLU_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGDLU_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGDLU_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGDLU_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGDLU_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGDLU_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGDLU_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGDLU_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGDLU_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGDLU_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGDLU_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGDLU_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGDLU_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGDLU_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGDLU_D1_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGDLU_D1_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGDLU_D1_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGDLU_D1_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGDLU_D1_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGDLU_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGDLU_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGDLU_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGDLU_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGDLU_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGDLU_D1_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGDLU_D1_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGDLU_D1_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGDLU_D1_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGDLU_D1_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGDLU_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGDLU_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGDLU_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGDLU_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGDLU_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGDLU_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGDLU_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGDLU_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGDLU_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGDLU_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
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

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGDLU_D2_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGDLU_D2_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGDLU_D2_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGDLU_D2_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGDLU_D2_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGDLU_D2_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGDLU_D2_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGDLU_D2_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGDLU_D2_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGDLU_D2_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGDLU_D2_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGDLU_D2_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGDLU_D2_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGDLU_D2_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGDLU_D2_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGDLU_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGDLU_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGDLU_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGDLU_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGDLU_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGDLU_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGDLU_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGDLU_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGDLU_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGDLU_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGDLU_D3_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGDLU_D3_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGDLU_D3_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGDLU_D3_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGDLU_D3_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGDLU_D3_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGDLU_D3_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGDLU_D3_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGDLU_D3_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGDLU_D3_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGDLU_D3_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGDLU_D3_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGDLU_D3_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGDLU_D3_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGDLU_D3_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGDLU_D3_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGDLU_D3_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGDLU_D3_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGDLU_D3_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGDLU_D3_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGDLU_D3_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGDLU_D3_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGDLU_D3_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGDLU_D3_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGDLU_D3_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D3_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1
    module procedure getUnifRandRNGFDD_D0_LK
        use pm_kind, only: LKC => LK
#include "pm_distUnif@routines.inc.F90"
    end procedure
#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGFLU_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGFLU_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGFLU_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGFLU_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGFLU_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGFLU_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGFLU_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGFLU_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGFLU_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGFLU_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGFLU_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGFLU_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGFLU_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGFLU_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGFLU_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGFLU_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGFLU_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGFLU_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGFLU_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGFLU_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGFLU_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGFLU_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGFLU_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGFLU_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGFLU_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGFLU_D1_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGFLU_D1_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGFLU_D1_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGFLU_D1_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGFLU_D1_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGFLU_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGFLU_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGFLU_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGFLU_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGFLU_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGFLU_D1_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGFLU_D1_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGFLU_D1_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGFLU_D1_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGFLU_D1_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGFLU_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGFLU_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGFLU_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGFLU_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGFLU_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGFLU_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGFLU_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGFLU_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGFLU_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGFLU_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
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

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGFLU_D2_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGFLU_D2_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGFLU_D2_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGFLU_D2_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGFLU_D2_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGFLU_D2_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGFLU_D2_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGFLU_D2_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGFLU_D2_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGFLU_D2_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGFLU_D2_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGFLU_D2_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGFLU_D2_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGFLU_D2_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGFLU_D2_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGFLU_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGFLU_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGFLU_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGFLU_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGFLU_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGFLU_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGFLU_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGFLU_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGFLU_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGFLU_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGFLU_D3_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGFLU_D3_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGFLU_D3_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGFLU_D3_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGFLU_D3_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGFLU_D3_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGFLU_D3_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGFLU_D3_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGFLU_D3_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGFLU_D3_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGFLU_D3_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGFLU_D3_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGFLU_D3_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGFLU_D3_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGFLU_D3_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGFLU_D3_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGFLU_D3_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGFLU_D3_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGFLU_D3_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGFLU_D3_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGFLU_D3_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGFLU_D3_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGFLU_D3_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGFLU_D3_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGFLU_D3_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D3_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SM64_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1
    module procedure getUnifRandRNGSDD_D0_LK
        use pm_kind, only: LKC => LK
#include "pm_distUnif@routines.inc.F90"
    end procedure
#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGSLU_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGSLU_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGSLU_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGSLU_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGSLU_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGSLU_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGSLU_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGSLU_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGSLU_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGSLU_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGSLU_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGSLU_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGSLU_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGSLU_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGSLU_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGSLU_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGSLU_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGSLU_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGSLU_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGSLU_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGSLU_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGSLU_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGSLU_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGSLU_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGSLU_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGSLU_D1_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGSLU_D1_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGSLU_D1_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGSLU_D1_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGSLU_D1_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGSLU_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGSLU_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGSLU_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGSLU_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGSLU_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGSLU_D1_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGSLU_D1_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGSLU_D1_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGSLU_D1_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGSLU_D1_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGSLU_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGSLU_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGSLU_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGSLU_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGSLU_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGSLU_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGSLU_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGSLU_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGSLU_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGSLU_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
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

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGSLU_D2_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGSLU_D2_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGSLU_D2_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGSLU_D2_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGSLU_D2_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGSLU_D2_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGSLU_D2_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGSLU_D2_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGSLU_D2_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGSLU_D2_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGSLU_D2_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGSLU_D2_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGSLU_D2_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGSLU_D2_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGSLU_D2_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGSLU_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGSLU_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGSLU_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGSLU_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGSLU_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGSLU_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGSLU_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGSLU_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGSLU_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGSLU_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGSLU_D3_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGSLU_D3_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGSLU_D3_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGSLU_D3_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGSLU_D3_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGSLU_D3_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGSLU_D3_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGSLU_D3_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGSLU_D3_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGSLU_D3_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGSLU_D3_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGSLU_D3_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGSLU_D3_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGSLU_D3_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGSLU_D3_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGSLU_D3_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGSLU_D3_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGSLU_D3_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGSLU_D3_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGSLU_D3_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGSLU_D3_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGSLU_D3_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGSLU_D3_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGSLU_D3_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGSLU_D3_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D3_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SM64_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define X256SSG_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1
    module procedure getUnifRandRNGGDD_D0_LK
        use pm_kind, only: LKC => LK
#include "pm_distUnif@routines.inc.F90"
    end procedure
#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGGLU_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGGLU_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGGLU_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGGLU_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGGLU_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGGLU_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGGLU_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGGLU_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGGLU_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGGLU_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGGLU_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGGLU_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGGLU_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGGLU_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGGLU_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGGLU_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGGLU_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGGLU_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGGLU_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGGLU_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGGLU_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGGLU_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGGLU_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGGLU_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGGLU_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGGLU_D1_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGGLU_D1_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGGLU_D1_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGGLU_D1_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGGLU_D1_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGGLU_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGGLU_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGGLU_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGGLU_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGGLU_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGGLU_D1_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGGLU_D1_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGGLU_D1_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGGLU_D1_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGGLU_D1_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGGLU_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGGLU_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGGLU_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGGLU_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGGLU_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGGLU_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGGLU_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGGLU_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGGLU_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGGLU_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
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

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGGLU_D2_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGGLU_D2_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGGLU_D2_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGGLU_D2_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGGLU_D2_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGGLU_D2_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGGLU_D2_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGGLU_D2_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGGLU_D2_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGGLU_D2_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGGLU_D2_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGGLU_D2_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGGLU_D2_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGGLU_D2_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGGLU_D2_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGGLU_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGGLU_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGGLU_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGGLU_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGGLU_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGGLU_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGGLU_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGGLU_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGGLU_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGGLU_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGGLU_D3_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGGLU_D3_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGGLU_D3_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGGLU_D3_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGGLU_D3_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGGLU_D3_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGGLU_D3_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGGLU_D3_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGGLU_D3_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGGLU_D3_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGGLU_D3_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGGLU_D3_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGGLU_D3_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGGLU_D3_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGGLU_D3_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGGLU_D3_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGGLU_D3_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGGLU_D3_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGGLU_D3_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGGLU_D3_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGGLU_D3_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGGLU_D3_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGGLU_D3_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGGLU_D3_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGGLU_D3_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D3_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef X256SSG_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define X256SSW_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1
    module procedure getUnifRandRNGXDD_D0_LK
        use pm_kind, only: LKC => LK
#include "pm_distUnif@routines.inc.F90"
    end procedure
#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGXLU_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGXLU_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGXLU_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGXLU_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGXLU_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGXLU_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGXLU_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGXLU_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGXLU_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGXLU_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGXLU_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGXLU_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGXLU_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGXLU_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGXLU_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGXLU_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGXLU_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGXLU_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGXLU_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGXLU_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGXLU_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGXLU_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGXLU_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGXLU_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGXLU_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGXLU_D1_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGXLU_D1_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGXLU_D1_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGXLU_D1_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGXLU_D1_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGXLU_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGXLU_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGXLU_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGXLU_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGXLU_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGXLU_D1_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGXLU_D1_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGXLU_D1_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGXLU_D1_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGXLU_D1_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGXLU_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGXLU_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGXLU_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGXLU_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGXLU_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGXLU_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGXLU_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGXLU_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGXLU_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGXLU_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
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

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGXLU_D2_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGXLU_D2_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGXLU_D2_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGXLU_D2_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGXLU_D2_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGXLU_D2_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGXLU_D2_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGXLU_D2_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGXLU_D2_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGXLU_D2_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGXLU_D2_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGXLU_D2_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGXLU_D2_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGXLU_D2_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGXLU_D2_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGXLU_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGXLU_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGXLU_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGXLU_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGXLU_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGXLU_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGXLU_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGXLU_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGXLU_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGXLU_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getUnifRandRNGXLU_D3_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getUnifRandRNGXLU_D3_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getUnifRandRNGXLU_D3_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getUnifRandRNGXLU_D3_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getUnifRandRNGXLU_D3_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getUnifRandRNGXLU_D3_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getUnifRandRNGXLU_D3_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getUnifRandRNGXLU_D3_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getUnifRandRNGXLU_D3_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getUnifRandRNGXLU_D3_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getUnifRandRNGXLU_D3_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getUnifRandRNGXLU_D3_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getUnifRandRNGXLU_D3_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getUnifRandRNGXLU_D3_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getUnifRandRNGXLU_D3_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getUnifRandRNGXLU_D3_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getUnifRandRNGXLU_D3_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getUnifRandRNGXLU_D3_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getUnifRandRNGXLU_D3_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getUnifRandRNGXLU_D3_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRandRNGXLU_D3_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRandRNGXLU_D3_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRandRNGXLU_D3_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRandRNGXLU_D3_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRandRNGXLU_D3_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D3_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef X256SSW_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getUnifRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setUnifRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1
#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGDDD_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGDDD_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGDDD_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGDDD_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGDDD_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGDDD_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGDDD_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGDDD_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGDDD_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGDDD_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGDDD_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGDDD_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGDDD_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGDDD_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGDDD_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGDDD_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGDDD_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGDDD_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGDDD_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGDDD_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGDDD_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGDDD_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGDDD_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGDDD_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGDDD_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGDLU_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGDLU_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGDLU_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGDLU_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGDLU_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGDLU_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGDLU_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGDLU_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGDLU_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGDLU_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGDLU_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGDLU_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGDLU_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGDLU_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGDLU_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGDLU_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGDLU_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGDLU_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGDLU_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGDLU_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGDLU_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGDLU_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGDLU_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGDLU_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGDLU_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGD_ENABLED
#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGFDD_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGFDD_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGFDD_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGFDD_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGFDD_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGFDD_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGFDD_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGFDD_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGFDD_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGFDD_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGFDD_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGFDD_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGFDD_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGFDD_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGFDD_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGFDD_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGFDD_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGFDD_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGFDD_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGFDD_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGFDD_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGFDD_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGFDD_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGFDD_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGFDD_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGFLU_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGFLU_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGFLU_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGFLU_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGFLU_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGFLU_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGFLU_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGFLU_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGFLU_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGFLU_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGFLU_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGFLU_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGFLU_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGFLU_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGFLU_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGFLU_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGFLU_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGFLU_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGFLU_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGFLU_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGFLU_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGFLU_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGFLU_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGFLU_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGFLU_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SM64_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGSDD_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGSDD_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGSDD_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGSDD_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGSDD_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGSDD_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGSDD_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGSDD_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGSDD_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGSDD_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGSDD_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGSDD_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGSDD_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGSDD_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGSDD_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGSDD_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGSDD_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGSDD_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGSDD_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGSDD_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGSDD_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGSDD_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGSDD_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGSDD_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGSDD_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGSDD_D1_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGSDD_D1_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGSDD_D1_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGSDD_D1_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGSDD_D1_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGSDD_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGSDD_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGSDD_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGSDD_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGSDD_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGSDD_D1_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGSDD_D1_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGSDD_D1_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGSDD_D1_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGSDD_D1_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGSDD_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGSDD_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGSDD_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGSDD_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGSDD_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGSDD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGSDD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGSDD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGSDD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGSDD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
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

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGSDD_D2_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGSDD_D2_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGSDD_D2_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGSDD_D2_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGSDD_D2_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGSDD_D2_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGSDD_D2_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGSDD_D2_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGSDD_D2_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGSDD_D2_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGSDD_D2_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGSDD_D2_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGSDD_D2_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGSDD_D2_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGSDD_D2_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGSDD_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGSDD_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGSDD_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGSDD_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGSDD_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGSDD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGSDD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGSDD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGSDD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGSDD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGSDD_D3_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGSDD_D3_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGSDD_D3_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGSDD_D3_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGSDD_D3_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGSDD_D3_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGSDD_D3_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGSDD_D3_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGSDD_D3_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGSDD_D3_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGSDD_D3_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGSDD_D3_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGSDD_D3_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGSDD_D3_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGSDD_D3_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGSDD_D3_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGSDD_D3_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGSDD_D3_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGSDD_D3_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGSDD_D3_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGSDD_D3_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGSDD_D3_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGSDD_D3_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGSDD_D3_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGSDD_D3_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D3_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGSLU_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGSLU_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGSLU_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGSLU_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGSLU_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGSLU_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGSLU_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGSLU_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGSLU_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGSLU_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGSLU_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGSLU_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGSLU_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGSLU_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGSLU_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGSLU_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGSLU_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGSLU_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGSLU_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGSLU_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGSLU_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGSLU_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGSLU_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGSLU_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGSLU_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGSLU_D1_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGSLU_D1_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGSLU_D1_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGSLU_D1_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGSLU_D1_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGSLU_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGSLU_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGSLU_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGSLU_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGSLU_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGSLU_D1_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGSLU_D1_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGSLU_D1_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGSLU_D1_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGSLU_D1_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGSLU_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGSLU_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGSLU_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGSLU_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGSLU_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGSLU_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGSLU_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGSLU_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGSLU_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGSLU_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
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

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGSLU_D2_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGSLU_D2_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGSLU_D2_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGSLU_D2_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGSLU_D2_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGSLU_D2_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGSLU_D2_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGSLU_D2_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGSLU_D2_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGSLU_D2_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGSLU_D2_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGSLU_D2_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGSLU_D2_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGSLU_D2_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGSLU_D2_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGSLU_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGSLU_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGSLU_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGSLU_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGSLU_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGSLU_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGSLU_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGSLU_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGSLU_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGSLU_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGSLU_D3_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGSLU_D3_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGSLU_D3_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGSLU_D3_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGSLU_D3_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGSLU_D3_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGSLU_D3_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGSLU_D3_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGSLU_D3_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGSLU_D3_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGSLU_D3_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGSLU_D3_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGSLU_D3_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGSLU_D3_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGSLU_D3_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGSLU_D3_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGSLU_D3_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGSLU_D3_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGSLU_D3_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGSLU_D3_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGSLU_D3_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGSLU_D3_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGSLU_D3_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGSLU_D3_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGSLU_D3_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D3_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SM64_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define X256SSG_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGGDD_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGGDD_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGGDD_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGGDD_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGGDD_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGGDD_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGGDD_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGGDD_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGGDD_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGGDD_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGGDD_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGGDD_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGGDD_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGGDD_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGGDD_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGGDD_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGGDD_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGGDD_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGGDD_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGGDD_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGGDD_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGGDD_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGGDD_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGGDD_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGGDD_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGGDD_D1_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGGDD_D1_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGGDD_D1_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGGDD_D1_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGGDD_D1_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGGDD_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGGDD_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGGDD_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGGDD_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGGDD_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGGDD_D1_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGGDD_D1_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGGDD_D1_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGGDD_D1_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGGDD_D1_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGGDD_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGGDD_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGGDD_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGGDD_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGGDD_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGGDD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGGDD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGGDD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGGDD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGGDD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
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

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGGDD_D2_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGGDD_D2_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGGDD_D2_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGGDD_D2_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGGDD_D2_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGGDD_D2_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGGDD_D2_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGGDD_D2_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGGDD_D2_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGGDD_D2_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGGDD_D2_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGGDD_D2_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGGDD_D2_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGGDD_D2_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGGDD_D2_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGGDD_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGGDD_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGGDD_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGGDD_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGGDD_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGGDD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGGDD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGGDD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGGDD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGGDD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGGDD_D3_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGGDD_D3_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGGDD_D3_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGGDD_D3_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGGDD_D3_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGGDD_D3_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGGDD_D3_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGGDD_D3_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGGDD_D3_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGGDD_D3_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGGDD_D3_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGGDD_D3_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGGDD_D3_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGGDD_D3_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGGDD_D3_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGGDD_D3_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGGDD_D3_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGGDD_D3_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGGDD_D3_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGGDD_D3_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGGDD_D3_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGGDD_D3_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGGDD_D3_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGGDD_D3_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGGDD_D3_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D3_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGGLU_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGGLU_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGGLU_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGGLU_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGGLU_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGGLU_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGGLU_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGGLU_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGGLU_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGGLU_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGGLU_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGGLU_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGGLU_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGGLU_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGGLU_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGGLU_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGGLU_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGGLU_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGGLU_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGGLU_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGGLU_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGGLU_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGGLU_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGGLU_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGGLU_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGGLU_D1_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGGLU_D1_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGGLU_D1_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGGLU_D1_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGGLU_D1_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGGLU_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGGLU_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGGLU_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGGLU_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGGLU_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGGLU_D1_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGGLU_D1_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGGLU_D1_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGGLU_D1_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGGLU_D1_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGGLU_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGGLU_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGGLU_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGGLU_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGGLU_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGGLU_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGGLU_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGGLU_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGGLU_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGGLU_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
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

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGGLU_D2_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGGLU_D2_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGGLU_D2_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGGLU_D2_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGGLU_D2_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGGLU_D2_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGGLU_D2_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGGLU_D2_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGGLU_D2_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGGLU_D2_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGGLU_D2_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGGLU_D2_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGGLU_D2_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGGLU_D2_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGGLU_D2_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGGLU_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGGLU_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGGLU_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGGLU_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGGLU_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGGLU_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGGLU_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGGLU_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGGLU_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGGLU_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGGLU_D3_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGGLU_D3_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGGLU_D3_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGGLU_D3_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGGLU_D3_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGGLU_D3_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGGLU_D3_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGGLU_D3_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGGLU_D3_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGGLU_D3_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGGLU_D3_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGGLU_D3_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGGLU_D3_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGGLU_D3_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGGLU_D3_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGGLU_D3_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGGLU_D3_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGGLU_D3_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGGLU_D3_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGGLU_D3_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGGLU_D3_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGGLU_D3_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGGLU_D3_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGGLU_D3_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGGLU_D3_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D3_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef X256SSG_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define X256SSW_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGXDD_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGXDD_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGXDD_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGXDD_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGXDD_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGXDD_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGXDD_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGXDD_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGXDD_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGXDD_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGXDD_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGXDD_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGXDD_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGXDD_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGXDD_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGXDD_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGXDD_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGXDD_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGXDD_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGXDD_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGXDD_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGXDD_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGXDD_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGXDD_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGXDD_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGXDD_D1_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGXDD_D1_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGXDD_D1_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGXDD_D1_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGXDD_D1_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGXDD_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGXDD_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGXDD_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGXDD_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGXDD_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGXDD_D1_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGXDD_D1_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGXDD_D1_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGXDD_D1_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGXDD_D1_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGXDD_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGXDD_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGXDD_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGXDD_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGXDD_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGXDD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGXDD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGXDD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGXDD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGXDD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
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

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGXDD_D2_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGXDD_D2_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGXDD_D2_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGXDD_D2_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGXDD_D2_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGXDD_D2_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGXDD_D2_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGXDD_D2_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGXDD_D2_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGXDD_D2_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGXDD_D2_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGXDD_D2_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGXDD_D2_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGXDD_D2_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGXDD_D2_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGXDD_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGXDD_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGXDD_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGXDD_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGXDD_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGXDD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGXDD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGXDD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGXDD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGXDD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGXDD_D3_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGXDD_D3_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGXDD_D3_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGXDD_D3_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGXDD_D3_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGXDD_D3_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGXDD_D3_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGXDD_D3_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGXDD_D3_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGXDD_D3_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGXDD_D3_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGXDD_D3_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGXDD_D3_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGXDD_D3_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGXDD_D3_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGXDD_D3_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGXDD_D3_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGXDD_D3_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGXDD_D3_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGXDD_D3_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGXDD_D3_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGXDD_D3_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGXDD_D3_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGXDD_D3_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGXDD_D3_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D3_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGXLU_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGXLU_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGXLU_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGXLU_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGXLU_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGXLU_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGXLU_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGXLU_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGXLU_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGXLU_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGXLU_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGXLU_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGXLU_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGXLU_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGXLU_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGXLU_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGXLU_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGXLU_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGXLU_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGXLU_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGXLU_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGXLU_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGXLU_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGXLU_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGXLU_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGXLU_D1_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGXLU_D1_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGXLU_D1_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGXLU_D1_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGXLU_D1_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGXLU_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGXLU_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGXLU_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGXLU_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGXLU_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGXLU_D1_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGXLU_D1_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGXLU_D1_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGXLU_D1_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGXLU_D1_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGXLU_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGXLU_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGXLU_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGXLU_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGXLU_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGXLU_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGXLU_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGXLU_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGXLU_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGXLU_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
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

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGXLU_D2_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGXLU_D2_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGXLU_D2_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGXLU_D2_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGXLU_D2_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGXLU_D2_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGXLU_D2_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGXLU_D2_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGXLU_D2_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGXLU_D2_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGXLU_D2_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGXLU_D2_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGXLU_D2_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGXLU_D2_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGXLU_D2_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGXLU_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGXLU_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGXLU_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGXLU_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGXLU_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGXLU_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGXLU_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGXLU_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGXLU_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGXLU_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setUnifRandRNGXLU_D3_SK5
        use pm_kind, only: SKC => SK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setUnifRandRNGXLU_D3_SK4
        use pm_kind, only: SKC => SK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setUnifRandRNGXLU_D3_SK3
        use pm_kind, only: SKC => SK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setUnifRandRNGXLU_D3_SK2
        use pm_kind, only: SKC => SK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setUnifRandRNGXLU_D3_SK1
        use pm_kind, only: SKC => SK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setUnifRandRNGXLU_D3_IK5
        use pm_kind, only: IKC => IK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setUnifRandRNGXLU_D3_IK4
        use pm_kind, only: IKC => IK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setUnifRandRNGXLU_D3_IK3
        use pm_kind, only: IKC => IK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setUnifRandRNGXLU_D3_IK2
        use pm_kind, only: IKC => IK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setUnifRandRNGXLU_D3_IK1
        use pm_kind, only: IKC => IK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setUnifRandRNGXLU_D3_LK5
        use pm_kind, only: LKC => LK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setUnifRandRNGXLU_D3_LK4
        use pm_kind, only: LKC => LK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setUnifRandRNGXLU_D3_LK3
        use pm_kind, only: LKC => LK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setUnifRandRNGXLU_D3_LK2
        use pm_kind, only: LKC => LK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setUnifRandRNGXLU_D3_LK1
        use pm_kind, only: LKC => LK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setUnifRandRNGXLU_D3_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setUnifRandRNGXLU_D3_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setUnifRandRNGXLU_D3_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setUnifRandRNGXLU_D3_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setUnifRandRNGXLU_D3_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRandRNGXLU_D3_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRandRNGXLU_D3_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRandRNGXLU_D3_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRandRNGXLU_D3_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRandRNGXLU_D3_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D3_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef X256SSW_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setUnifRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines