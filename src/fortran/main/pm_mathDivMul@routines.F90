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
!>  This file contains procedure implementations of [pm_mathDivMul](@ref pm_mathDivMul).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_mathDivMul) routines ! LCOV_EXCL_LINE

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

#define getDivMul_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unary_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getDivMulUnary_IK5
        use pm_kind, only: IKG => IK5
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getDivMulUnary_IK4
        use pm_kind, only: IKG => IK4
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getDivMulUnary_IK3
        use pm_kind, only: IKG => IK3
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getDivMulUnary_IK2
        use pm_kind, only: IKG => IK2
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getDivMulUnary_IK1
        use pm_kind, only: IKG => IK1
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDivMulUnary_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDivMulUnary_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDivMulUnary_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDivMulUnary_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDivMulUnary_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getDivMulUnary_CK5
        use pm_kind, only: CKG => CK5
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getDivMulUnary_CK4
        use pm_kind, only: CKG => CK4
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getDivMulUnary_CK3
        use pm_kind, only: CKG => CK3
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getDivMulUnary_CK2
        use pm_kind, only: CKG => CK2
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getDivMulUnary_CK1
        use pm_kind, only: CKG => CK1
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Unary_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Binary_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getDivMulBinary_IK5
        use pm_kind, only: IKG => IK5
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getDivMulBinary_IK4
        use pm_kind, only: IKG => IK4
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getDivMulBinary_IK3
        use pm_kind, only: IKG => IK3
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getDivMulBinary_IK2
        use pm_kind, only: IKG => IK2
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getDivMulBinary_IK1
        use pm_kind, only: IKG => IK1
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDivMulBinary_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDivMulBinary_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDivMulBinary_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDivMulBinary_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDivMulBinary_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getDivMulBinary_CK5
        use pm_kind, only: CKG => CK5
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getDivMulBinary_CK4
        use pm_kind, only: CKG => CK4
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getDivMulBinary_CK3
        use pm_kind, only: CKG => CK3
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getDivMulBinary_CK2
        use pm_kind, only: CKG => CK2
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getDivMulBinary_CK1
        use pm_kind, only: CKG => CK1
#include "pm_mathDivMul@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Binary_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDivMul_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
