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
!>  This file contains procedure implementations of [pm_except](@ref pm_except).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_except) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_negative_inf
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_negative

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isAddOutflow_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure isAddOutflow_IK5
        use pm_kind, only: IKC => IK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure isAddOutflow_IK4
        use pm_kind, only: IKC => IK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure isAddOutflow_IK3
        use pm_kind, only: IKC => IK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure isAddOutflow_IK2
        use pm_kind, only: IKC => IK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure isAddOutflow_IK1
        use pm_kind, only: IKC => IK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isAddOutflow_CK5
        use pm_kind, only: CKC => CK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isAddOutflow_CK4
        use pm_kind, only: CKC => CK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isAddOutflow_CK3
        use pm_kind, only: CKC => CK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isAddOutflow_CK2
        use pm_kind, only: CKC => CK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isAddOutflow_CK1
        use pm_kind, only: CKC => CK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isAddOutflow_RK5
        use pm_kind, only: RKC => RK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isAddOutflow_RK4
        use pm_kind, only: RKC => RK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isAddOutflow_RK3
        use pm_kind, only: RKC => RK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isAddOutflow_RK2
        use pm_kind, only: RKC => RK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isAddOutflow_RK1
        use pm_kind, only: RKC => RK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isAddOutflow_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isAddOutflowNeg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure isAddOutflowNeg_IK5
        use pm_kind, only: IKC => IK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure isAddOutflowNeg_IK4
        use pm_kind, only: IKC => IK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure isAddOutflowNeg_IK3
        use pm_kind, only: IKC => IK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure isAddOutflowNeg_IK2
        use pm_kind, only: IKC => IK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure isAddOutflowNeg_IK1
        use pm_kind, only: IKC => IK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isAddOutflowNeg_CK5
        use pm_kind, only: CKC => CK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isAddOutflowNeg_CK4
        use pm_kind, only: CKC => CK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isAddOutflowNeg_CK3
        use pm_kind, only: CKC => CK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isAddOutflowNeg_CK2
        use pm_kind, only: CKC => CK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isAddOutflowNeg_CK1
        use pm_kind, only: CKC => CK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isAddOutflowNeg_RK5
        use pm_kind, only: RKC => RK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isAddOutflowNeg_RK4
        use pm_kind, only: RKC => RK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isAddOutflowNeg_RK3
        use pm_kind, only: RKC => RK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isAddOutflowNeg_RK2
        use pm_kind, only: RKC => RK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isAddOutflowNeg_RK1
        use pm_kind, only: RKC => RK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isAddOutflowNeg_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isAddOutflowPos_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure isAddOutflowPos_IK5
        use pm_kind, only: IKC => IK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure isAddOutflowPos_IK4
        use pm_kind, only: IKC => IK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure isAddOutflowPos_IK3
        use pm_kind, only: IKC => IK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure isAddOutflowPos_IK2
        use pm_kind, only: IKC => IK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure isAddOutflowPos_IK1
        use pm_kind, only: IKC => IK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isAddOutflowPos_CK5
        use pm_kind, only: CKC => CK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isAddOutflowPos_CK4
        use pm_kind, only: CKC => CK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isAddOutflowPos_CK3
        use pm_kind, only: CKC => CK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isAddOutflowPos_CK2
        use pm_kind, only: CKC => CK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isAddOutflowPos_CK1
        use pm_kind, only: CKC => CK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isAddOutflowPos_RK5
        use pm_kind, only: RKC => RK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isAddOutflowPos_RK4
        use pm_kind, only: RKC => RK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isAddOutflowPos_RK3
        use pm_kind, only: RKC => RK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isAddOutflowPos_RK2
        use pm_kind, only: RKC => RK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isAddOutflowPos_RK1
        use pm_kind, only: RKC => RK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isAddOutflowPos_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isInf_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isInf_RK5
        use pm_kind, only: RKC => RK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isInf_RK4
        use pm_kind, only: RKC => RK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isInf_RK3
        use pm_kind, only: RKC => RK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isInf_RK2
        use pm_kind, only: RKC => RK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isInf_RK1
        use pm_kind, only: RKC => RK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isInf_CK5
        use pm_kind, only: CKC => CK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isInf_CK4
        use pm_kind, only: CKC => CK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isInf_CK3
        use pm_kind, only: CKC => CK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isInf_CK2
        use pm_kind, only: CKC => CK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isInf_CK1
        use pm_kind, only: CKC => CK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isInf_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isInfPos_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isInfPos_RK5
        use pm_kind, only: RKC => RK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isInfPos_RK4
        use pm_kind, only: RKC => RK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isInfPos_RK3
        use pm_kind, only: RKC => RK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isInfPos_RK2
        use pm_kind, only: RKC => RK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isInfPos_RK1
        use pm_kind, only: RKC => RK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isInfPos_CK5
        use pm_kind, only: CKC => CK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isInfPos_CK4
        use pm_kind, only: CKC => CK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isInfPos_CK3
        use pm_kind, only: CKC => CK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isInfPos_CK2
        use pm_kind, only: CKC => CK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isInfPos_CK1
        use pm_kind, only: CKC => CK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isInfPos_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getInfPos_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getInfPos_RK5
        use pm_kind, only: RKC => RK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getInfPos_RK4
        use pm_kind, only: RKC => RK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getInfPos_RK3
        use pm_kind, only: RKC => RK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getInfPos_RK2
        use pm_kind, only: RKC => RK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getInfPos_RK1
        use pm_kind, only: RKC => RK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getInfPos_CK5
        use pm_kind, only: CKC => CK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getInfPos_CK4
        use pm_kind, only: CKC => CK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getInfPos_CK3
        use pm_kind, only: CKC => CK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getInfPos_CK2
        use pm_kind, only: CKC => CK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getInfPos_CK1
        use pm_kind, only: CKC => CK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getInfPos_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setInfPos_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setInfPos_RK5
        use pm_kind, only: RKC => RK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setInfPos_RK4
        use pm_kind, only: RKC => RK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setInfPos_RK3
        use pm_kind, only: RKC => RK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setInfPos_RK2
        use pm_kind, only: RKC => RK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setInfPos_RK1
        use pm_kind, only: RKC => RK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setInfPos_CK5
        use pm_kind, only: CKC => CK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setInfPos_CK4
        use pm_kind, only: CKC => CK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setInfPos_CK3
        use pm_kind, only: CKC => CK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setInfPos_CK2
        use pm_kind, only: CKC => CK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setInfPos_CK1
        use pm_kind, only: CKC => CK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setInfPos_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isInfNeg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isInfNeg_RK5
        use pm_kind, only: RKC => RK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isInfNeg_RK4
        use pm_kind, only: RKC => RK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isInfNeg_RK3
        use pm_kind, only: RKC => RK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isInfNeg_RK2
        use pm_kind, only: RKC => RK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isInfNeg_RK1
        use pm_kind, only: RKC => RK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isInfNeg_CK5
        use pm_kind, only: CKC => CK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isInfNeg_CK4
        use pm_kind, only: CKC => CK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isInfNeg_CK3
        use pm_kind, only: CKC => CK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isInfNeg_CK2
        use pm_kind, only: CKC => CK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isInfNeg_CK1
        use pm_kind, only: CKC => CK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isInfNeg_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getInfNeg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getInfNeg_RK5
        use pm_kind, only: RKC => RK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getInfNeg_RK4
        use pm_kind, only: RKC => RK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getInfNeg_RK3
        use pm_kind, only: RKC => RK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getInfNeg_RK2
        use pm_kind, only: RKC => RK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getInfNeg_RK1
        use pm_kind, only: RKC => RK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getInfNeg_CK5
        use pm_kind, only: CKC => CK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getInfNeg_CK4
        use pm_kind, only: CKC => CK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getInfNeg_CK3
        use pm_kind, only: CKC => CK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getInfNeg_CK2
        use pm_kind, only: CKC => CK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getInfNeg_CK1
        use pm_kind, only: CKC => CK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getInfNeg_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setInfNeg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setInfNeg_RK5
        use pm_kind, only: RKC => RK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setInfNeg_RK4
        use pm_kind, only: RKC => RK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setInfNeg_RK3
        use pm_kind, only: RKC => RK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setInfNeg_RK2
        use pm_kind, only: RKC => RK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setInfNeg_RK1
        use pm_kind, only: RKC => RK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setInfNeg_CK5
        use pm_kind, only: CKC => CK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setInfNeg_CK4
        use pm_kind, only: CKC => CK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setInfNeg_CK3
        use pm_kind, only: CKC => CK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setInfNeg_CK2
        use pm_kind, only: CKC => CK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setInfNeg_CK1
        use pm_kind, only: CKC => CK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setInfNeg_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isNAN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IEEE_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isNANIEEE_CK5
        use pm_kind, only: CKC => CK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isNANIEEE_CK4
        use pm_kind, only: CKC => CK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isNANIEEE_CK3
        use pm_kind, only: CKC => CK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isNANIEEE_CK2
        use pm_kind, only: CKC => CK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isNANIEEE_CK1
        use pm_kind, only: CKC => CK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isNANIEEE_RK5
        use pm_kind, only: RKC => RK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isNANIEEE_RK4
        use pm_kind, only: RKC => RK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isNANIEEE_RK3
        use pm_kind, only: RKC => RK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isNANIEEE_RK2
        use pm_kind, only: RKC => RK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isNANIEEE_RK1
        use pm_kind, only: RKC => RK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IEEE_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XNEQ_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isNANXNEQ_CK5
        use pm_kind, only: CKC => CK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isNANXNEQ_CK4
        use pm_kind, only: CKC => CK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isNANXNEQ_CK3
        use pm_kind, only: CKC => CK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isNANXNEQ_CK2
        use pm_kind, only: CKC => CK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isNANXNEQ_CK1
        use pm_kind, only: CKC => CK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isNANXNEQ_RK5
        use pm_kind, only: RKC => RK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isNANXNEQ_RK4
        use pm_kind, only: RKC => RK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isNANXNEQ_RK3
        use pm_kind, only: RKC => RK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isNANXNEQ_RK2
        use pm_kind, only: RKC => RK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isNANXNEQ_RK1
        use pm_kind, only: RKC => RK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XNEQ_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isNAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNAN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getNAN_CK5
        use pm_kind, only: CKC => CK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getNAN_CK4
        use pm_kind, only: CKC => CK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getNAN_CK3
        use pm_kind, only: CKC => CK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getNAN_CK2
        use pm_kind, only: CKC => CK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getNAN_CK1
        use pm_kind, only: CKC => CK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getNAN_RK5
        use pm_kind, only: RKC => RK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getNAN_RK4
        use pm_kind, only: RKC => RK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getNAN_RK3
        use pm_kind, only: RKC => RK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getNAN_RK2
        use pm_kind, only: RKC => RK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getNAN_RK1
        use pm_kind, only: RKC => RK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setNAN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setNAN_CK5
        use pm_kind, only: CKC => CK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setNAN_CK4
        use pm_kind, only: CKC => CK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setNAN_CK3
        use pm_kind, only: CKC => CK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setNAN_CK2
        use pm_kind, only: CKC => CK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setNAN_CK1
        use pm_kind, only: CKC => CK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNAN_RK5
        use pm_kind, only: RKC => RK5
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNAN_RK4
        use pm_kind, only: RKC => RK4
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNAN_RK3
        use pm_kind, only: RKC => RK3
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNAN_RK2
        use pm_kind, only: RKC => RK2
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNAN_RK1
        use pm_kind, only: RKC => RK1
#include "pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setNAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines
