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
!>  This file contains procedure implementations of [pm_matrixMulAdd](@ref pm_matrixMulAdd).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_matrixMulAdd) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    ! \bug Bypass Intel `ifort` 2022 compiler bug for too many use statements in submodule procedures.
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    !use pm_arrayReverse, only: getReversed
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

#define setMatMulAdd_ENABLED 1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define gemv_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SFA_ENABLED 1
#define SFB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ASS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TNB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TNA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemv_ASS_SFA_SFB_TNA_TNB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemv_ASS_SFA_SFB_TNA_TNB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemv_ASS_SFA_SFB_TNA_TNB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemv_ASS_SFA_SFB_TNA_TNB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemv_ASS_SFA_SFB_TNA_TNB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemv_ASS_SFA_SFB_TNA_TNB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemv_ASS_SFA_SFB_TNA_TNB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemv_ASS_SFA_SFB_TNA_TNB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemv_ASS_SFA_SFB_TNA_TNB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemv_ASS_SFA_SFB_TNA_TNB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemv_ASS_SFA_SFB_TNA_TNB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemv_ASS_SFA_SFB_TNA_TNB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemv_ASS_SFA_SFB_TNA_TNB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemv_ASS_SFA_SFB_TNA_TNB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemv_ASS_SFA_SFB_TNA_TNB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TNA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemv_ASS_SFA_SFB_TSA_TNB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemv_ASS_SFA_SFB_TSA_TNB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemv_ASS_SFA_SFB_TSA_TNB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemv_ASS_SFA_SFB_TSA_TNB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemv_ASS_SFA_SFB_TSA_TNB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemv_ASS_SFA_SFB_TSA_TNB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemv_ASS_SFA_SFB_TSA_TNB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemv_ASS_SFA_SFB_TSA_TNB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemv_ASS_SFA_SFB_TSA_TNB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemv_ASS_SFA_SFB_TSA_TNB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemv_ASS_SFA_SFB_TSA_TNB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemv_ASS_SFA_SFB_TSA_TNB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemv_ASS_SFA_SFB_TSA_TNB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemv_ASS_SFA_SFB_TSA_TNB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemv_ASS_SFA_SFB_TSA_TNB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define THA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemv_ASS_SFA_SFB_THA_TNB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemv_ASS_SFA_SFB_THA_TNB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemv_ASS_SFA_SFB_THA_TNB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemv_ASS_SFA_SFB_THA_TNB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemv_ASS_SFA_SFB_THA_TNB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemv_ASS_SFA_SFB_THA_TNB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemv_ASS_SFA_SFB_THA_TNB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemv_ASS_SFA_SFB_THA_TNB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemv_ASS_SFA_SFB_THA_TNB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemv_ASS_SFA_SFB_THA_TNB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemv_ASS_SFA_SFB_THA_TNB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemv_ASS_SFA_SFB_THA_TNB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemv_ASS_SFA_SFB_THA_TNB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemv_ASS_SFA_SFB_THA_TNB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemv_ASS_SFA_SFB_THA_TNB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef THA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TNB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ASS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define EXP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TNB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TNA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemv_EXP_SFA_SFB_TNA_TNB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemv_EXP_SFA_SFB_TNA_TNB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemv_EXP_SFA_SFB_TNA_TNB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemv_EXP_SFA_SFB_TNA_TNB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemv_EXP_SFA_SFB_TNA_TNB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemv_EXP_SFA_SFB_TNA_TNB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemv_EXP_SFA_SFB_TNA_TNB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemv_EXP_SFA_SFB_TNA_TNB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemv_EXP_SFA_SFB_TNA_TNB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemv_EXP_SFA_SFB_TNA_TNB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemv_EXP_SFA_SFB_TNA_TNB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemv_EXP_SFA_SFB_TNA_TNB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemv_EXP_SFA_SFB_TNA_TNB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemv_EXP_SFA_SFB_TNA_TNB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemv_EXP_SFA_SFB_TNA_TNB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TNA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemv_EXP_SFA_SFB_TSA_TNB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemv_EXP_SFA_SFB_TSA_TNB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemv_EXP_SFA_SFB_TSA_TNB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemv_EXP_SFA_SFB_TSA_TNB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemv_EXP_SFA_SFB_TSA_TNB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemv_EXP_SFA_SFB_TSA_TNB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemv_EXP_SFA_SFB_TSA_TNB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemv_EXP_SFA_SFB_TSA_TNB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemv_EXP_SFA_SFB_TSA_TNB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemv_EXP_SFA_SFB_TSA_TNB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemv_EXP_SFA_SFB_TSA_TNB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemv_EXP_SFA_SFB_TSA_TNB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemv_EXP_SFA_SFB_TSA_TNB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemv_EXP_SFA_SFB_TSA_TNB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemv_EXP_SFA_SFB_TSA_TNB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define THA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemv_EXP_SFA_SFB_THA_TNB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemv_EXP_SFA_SFB_THA_TNB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemv_EXP_SFA_SFB_THA_TNB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemv_EXP_SFA_SFB_THA_TNB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemv_EXP_SFA_SFB_THA_TNB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemv_EXP_SFA_SFB_THA_TNB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemv_EXP_SFA_SFB_THA_TNB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemv_EXP_SFA_SFB_THA_TNB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemv_EXP_SFA_SFB_THA_TNB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemv_EXP_SFA_SFB_THA_TNB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemv_EXP_SFA_SFB_THA_TNB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemv_EXP_SFA_SFB_THA_TNB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemv_EXP_SFA_SFB_THA_TNB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemv_EXP_SFA_SFB_THA_TNB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemv_EXP_SFA_SFB_THA_TNB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef THA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TNB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef EXP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SFA_ENABLED
#undef SFB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef gemv_ENABLED
#endif
!FORTRAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define spmv_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CSA_ENABLED 1
#define CNB_ENABLED 1
#define SFB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ASS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure spmv_ASS_CSA_SUA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure spmv_ASS_CSA_SUA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure spmv_ASS_CSA_SUA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure spmv_ASS_CSA_SUA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure spmv_ASS_CSA_SUA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure spmv_ASS_CSA_SUA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure spmv_ASS_CSA_SUA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure spmv_ASS_CSA_SUA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure spmv_ASS_CSA_SUA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure spmv_ASS_CSA_SUA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure spmv_ASS_CSA_SUA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure spmv_ASS_CSA_SUA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure spmv_ASS_CSA_SUA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure spmv_ASS_CSA_SUA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure spmv_ASS_CSA_SUA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure spmv_ASS_CSA_SLA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure spmv_ASS_CSA_SLA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure spmv_ASS_CSA_SLA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure spmv_ASS_CSA_SLA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure spmv_ASS_CSA_SLA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure spmv_ASS_CSA_SLA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure spmv_ASS_CSA_SLA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure spmv_ASS_CSA_SLA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure spmv_ASS_CSA_SLA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure spmv_ASS_CSA_SLA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure spmv_ASS_CSA_SLA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure spmv_ASS_CSA_SLA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure spmv_ASS_CSA_SLA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure spmv_ASS_CSA_SLA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure spmv_ASS_CSA_SLA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ASS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define EXP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure spmv_EXP_CSA_SUA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure spmv_EXP_CSA_SUA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure spmv_EXP_CSA_SUA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure spmv_EXP_CSA_SUA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure spmv_EXP_CSA_SUA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure spmv_EXP_CSA_SUA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure spmv_EXP_CSA_SUA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure spmv_EXP_CSA_SUA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure spmv_EXP_CSA_SUA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure spmv_EXP_CSA_SUA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure spmv_EXP_CSA_SUA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure spmv_EXP_CSA_SUA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure spmv_EXP_CSA_SUA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure spmv_EXP_CSA_SUA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure spmv_EXP_CSA_SUA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure spmv_EXP_CSA_SLA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure spmv_EXP_CSA_SLA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure spmv_EXP_CSA_SLA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure spmv_EXP_CSA_SLA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure spmv_EXP_CSA_SLA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure spmv_EXP_CSA_SLA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure spmv_EXP_CSA_SLA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure spmv_EXP_CSA_SLA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure spmv_EXP_CSA_SLA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure spmv_EXP_CSA_SLA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure spmv_EXP_CSA_SLA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure spmv_EXP_CSA_SLA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure spmv_EXP_CSA_SLA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure spmv_EXP_CSA_SLA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure spmv_EXP_CSA_SLA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef EXP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CNB_ENABLED
#undef SFB_ENABLED
#undef CSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef spmv_ENABLED
#endif
!FORTRAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define hpmv_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CHA_ENABLED 1
#define CNB_ENABLED 1
#define SFB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ASS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure hpmv_ASS_CHA_SUA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure hpmv_ASS_CHA_SUA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure hpmv_ASS_CHA_SUA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure hpmv_ASS_CHA_SUA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure hpmv_ASS_CHA_SUA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure hpmv_ASS_CHA_SUA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure hpmv_ASS_CHA_SUA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure hpmv_ASS_CHA_SUA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure hpmv_ASS_CHA_SUA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure hpmv_ASS_CHA_SUA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure hpmv_ASS_CHA_SUA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure hpmv_ASS_CHA_SUA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure hpmv_ASS_CHA_SUA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure hpmv_ASS_CHA_SUA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure hpmv_ASS_CHA_SUA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure hpmv_ASS_CHA_SLA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure hpmv_ASS_CHA_SLA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure hpmv_ASS_CHA_SLA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure hpmv_ASS_CHA_SLA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure hpmv_ASS_CHA_SLA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure hpmv_ASS_CHA_SLA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure hpmv_ASS_CHA_SLA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure hpmv_ASS_CHA_SLA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure hpmv_ASS_CHA_SLA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure hpmv_ASS_CHA_SLA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure hpmv_ASS_CHA_SLA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure hpmv_ASS_CHA_SLA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure hpmv_ASS_CHA_SLA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure hpmv_ASS_CHA_SLA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure hpmv_ASS_CHA_SLA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ASS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define EXP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure hpmv_EXP_CHA_SUA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure hpmv_EXP_CHA_SUA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure hpmv_EXP_CHA_SUA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure hpmv_EXP_CHA_SUA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure hpmv_EXP_CHA_SUA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure hpmv_EXP_CHA_SUA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure hpmv_EXP_CHA_SUA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure hpmv_EXP_CHA_SUA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure hpmv_EXP_CHA_SUA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure hpmv_EXP_CHA_SUA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure hpmv_EXP_CHA_SUA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure hpmv_EXP_CHA_SUA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure hpmv_EXP_CHA_SUA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure hpmv_EXP_CHA_SUA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure hpmv_EXP_CHA_SUA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure hpmv_EXP_CHA_SLA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure hpmv_EXP_CHA_SLA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure hpmv_EXP_CHA_SLA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure hpmv_EXP_CHA_SLA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure hpmv_EXP_CHA_SLA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure hpmv_EXP_CHA_SLA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure hpmv_EXP_CHA_SLA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure hpmv_EXP_CHA_SLA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure hpmv_EXP_CHA_SLA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure hpmv_EXP_CHA_SLA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure hpmv_EXP_CHA_SLA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure hpmv_EXP_CHA_SLA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure hpmv_EXP_CHA_SLA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure hpmv_EXP_CHA_SLA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure hpmv_EXP_CHA_SLA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef EXP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CNB_ENABLED
#undef SFB_ENABLED
#undef CHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef hpmv_ENABLED
#endif
!FORTRAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define symv_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CSA_ENABLED 1
#define CNB_ENABLED 1
#define SFB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ASS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure symv_ASS_CSA_SUA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure symv_ASS_CSA_SUA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure symv_ASS_CSA_SUA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure symv_ASS_CSA_SUA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure symv_ASS_CSA_SUA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure symv_ASS_CSA_SUA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure symv_ASS_CSA_SUA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure symv_ASS_CSA_SUA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure symv_ASS_CSA_SUA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure symv_ASS_CSA_SUA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure symv_ASS_CSA_SUA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure symv_ASS_CSA_SUA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure symv_ASS_CSA_SUA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure symv_ASS_CSA_SUA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure symv_ASS_CSA_SUA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure symv_ASS_CSA_SLA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure symv_ASS_CSA_SLA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure symv_ASS_CSA_SLA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure symv_ASS_CSA_SLA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure symv_ASS_CSA_SLA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure symv_ASS_CSA_SLA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure symv_ASS_CSA_SLA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure symv_ASS_CSA_SLA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure symv_ASS_CSA_SLA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure symv_ASS_CSA_SLA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure symv_ASS_CSA_SLA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure symv_ASS_CSA_SLA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure symv_ASS_CSA_SLA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure symv_ASS_CSA_SLA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure symv_ASS_CSA_SLA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ASS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define EXP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure symv_EXP_CSA_SUA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure symv_EXP_CSA_SUA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure symv_EXP_CSA_SUA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure symv_EXP_CSA_SUA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure symv_EXP_CSA_SUA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure symv_EXP_CSA_SUA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure symv_EXP_CSA_SUA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure symv_EXP_CSA_SUA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure symv_EXP_CSA_SUA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure symv_EXP_CSA_SUA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure symv_EXP_CSA_SUA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure symv_EXP_CSA_SUA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure symv_EXP_CSA_SUA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure symv_EXP_CSA_SUA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure symv_EXP_CSA_SUA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure symv_EXP_CSA_SLA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure symv_EXP_CSA_SLA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure symv_EXP_CSA_SLA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure symv_EXP_CSA_SLA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure symv_EXP_CSA_SLA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure symv_EXP_CSA_SLA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure symv_EXP_CSA_SLA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure symv_EXP_CSA_SLA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure symv_EXP_CSA_SLA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure symv_EXP_CSA_SLA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure symv_EXP_CSA_SLA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure symv_EXP_CSA_SLA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure symv_EXP_CSA_SLA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure symv_EXP_CSA_SLA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure symv_EXP_CSA_SLA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef EXP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CNB_ENABLED
#undef SFB_ENABLED
#undef CSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef symv_ENABLED
#endif
!FORTRAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define hemv_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CHA_ENABLED 1
#define CNB_ENABLED 1
#define SFB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ASS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure hemv_ASS_CHA_SUA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure hemv_ASS_CHA_SUA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure hemv_ASS_CHA_SUA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure hemv_ASS_CHA_SUA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure hemv_ASS_CHA_SUA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure hemv_ASS_CHA_SUA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure hemv_ASS_CHA_SUA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure hemv_ASS_CHA_SUA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure hemv_ASS_CHA_SUA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure hemv_ASS_CHA_SUA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure hemv_ASS_CHA_SUA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure hemv_ASS_CHA_SUA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure hemv_ASS_CHA_SUA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure hemv_ASS_CHA_SUA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure hemv_ASS_CHA_SUA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure hemv_ASS_CHA_SLA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure hemv_ASS_CHA_SLA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure hemv_ASS_CHA_SLA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure hemv_ASS_CHA_SLA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure hemv_ASS_CHA_SLA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure hemv_ASS_CHA_SLA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure hemv_ASS_CHA_SLA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure hemv_ASS_CHA_SLA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure hemv_ASS_CHA_SLA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure hemv_ASS_CHA_SLA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure hemv_ASS_CHA_SLA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure hemv_ASS_CHA_SLA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure hemv_ASS_CHA_SLA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure hemv_ASS_CHA_SLA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure hemv_ASS_CHA_SLA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ASS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define EXP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure hemv_EXP_CHA_SUA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure hemv_EXP_CHA_SUA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure hemv_EXP_CHA_SUA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure hemv_EXP_CHA_SUA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure hemv_EXP_CHA_SUA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure hemv_EXP_CHA_SUA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure hemv_EXP_CHA_SUA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure hemv_EXP_CHA_SUA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure hemv_EXP_CHA_SUA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure hemv_EXP_CHA_SUA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure hemv_EXP_CHA_SUA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure hemv_EXP_CHA_SUA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure hemv_EXP_CHA_SUA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure hemv_EXP_CHA_SUA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure hemv_EXP_CHA_SUA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure hemv_EXP_CHA_SLA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure hemv_EXP_CHA_SLA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure hemv_EXP_CHA_SLA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure hemv_EXP_CHA_SLA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure hemv_EXP_CHA_SLA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure hemv_EXP_CHA_SLA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure hemv_EXP_CHA_SLA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure hemv_EXP_CHA_SLA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure hemv_EXP_CHA_SLA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure hemv_EXP_CHA_SLA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure hemv_EXP_CHA_SLA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure hemv_EXP_CHA_SLA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure hemv_EXP_CHA_SLA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure hemv_EXP_CHA_SLA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure hemv_EXP_CHA_SLA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef EXP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CNB_ENABLED
#undef SFB_ENABLED
#undef CHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef hemv_ENABLED
#endif
!FORTRAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define symm_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ASS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CNB_ENABLED 1
#define SFB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure symm_ASS_CSA_SUA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure symm_ASS_CSA_SUA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure symm_ASS_CSA_SUA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure symm_ASS_CSA_SUA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure symm_ASS_CSA_SUA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure symm_ASS_CSA_SUA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure symm_ASS_CSA_SUA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure symm_ASS_CSA_SUA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure symm_ASS_CSA_SUA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure symm_ASS_CSA_SUA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure symm_ASS_CSA_SUA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure symm_ASS_CSA_SUA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure symm_ASS_CSA_SUA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure symm_ASS_CSA_SUA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure symm_ASS_CSA_SUA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure symm_ASS_CSA_SLA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure symm_ASS_CSA_SLA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure symm_ASS_CSA_SLA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure symm_ASS_CSA_SLA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure symm_ASS_CSA_SLA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure symm_ASS_CSA_SLA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure symm_ASS_CSA_SLA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure symm_ASS_CSA_SLA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure symm_ASS_CSA_SLA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure symm_ASS_CSA_SLA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure symm_ASS_CSA_SLA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure symm_ASS_CSA_SLA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure symm_ASS_CSA_SLA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure symm_ASS_CSA_SLA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure symm_ASS_CSA_SLA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CNB_ENABLED
#undef SFB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CNA_ENABLED 1
#define SFA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SUB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SUB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SUB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SUB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SUB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SUB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SUB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SUB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SUB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SUB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SUB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SUB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SUB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SUB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SUB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SLB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SLB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SLB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SLB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SLB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SLB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SLB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SLB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SLB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SLB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SLB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SLB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SLB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SLB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure symm_ASS_CNA_SFA_CSB_SLB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CNA_ENABLED
#undef SFA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ASS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define EXP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CNB_ENABLED 1
#define SFB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure symm_EXP_CSA_SUA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure symm_EXP_CSA_SUA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure symm_EXP_CSA_SUA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure symm_EXP_CSA_SUA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure symm_EXP_CSA_SUA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure symm_EXP_CSA_SUA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure symm_EXP_CSA_SUA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure symm_EXP_CSA_SUA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure symm_EXP_CSA_SUA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure symm_EXP_CSA_SUA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure symm_EXP_CSA_SUA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure symm_EXP_CSA_SUA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure symm_EXP_CSA_SUA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure symm_EXP_CSA_SUA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure symm_EXP_CSA_SUA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure symm_EXP_CSA_SLA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure symm_EXP_CSA_SLA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure symm_EXP_CSA_SLA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure symm_EXP_CSA_SLA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure symm_EXP_CSA_SLA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure symm_EXP_CSA_SLA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure symm_EXP_CSA_SLA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure symm_EXP_CSA_SLA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure symm_EXP_CSA_SLA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure symm_EXP_CSA_SLA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure symm_EXP_CSA_SLA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure symm_EXP_CSA_SLA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure symm_EXP_CSA_SLA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure symm_EXP_CSA_SLA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure symm_EXP_CSA_SLA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CNB_ENABLED
#undef SFB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CNA_ENABLED 1
#define SFA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SUB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SUB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SUB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SUB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SUB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SUB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SUB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SUB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SUB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SUB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SUB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SUB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SUB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SUB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SUB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SLB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SLB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SLB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SLB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SLB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SLB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SLB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SLB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SLB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SLB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SLB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SLB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SLB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SLB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure symm_EXP_CNA_SFA_CSB_SLB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CNA_ENABLED
#undef SFA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef EXP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef symm_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define hemm_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ASS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CNB_ENABLED 1
#define SFB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure hemm_ASS_CHA_SUA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure hemm_ASS_CHA_SUA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure hemm_ASS_CHA_SUA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure hemm_ASS_CHA_SUA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure hemm_ASS_CHA_SUA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure hemm_ASS_CHA_SUA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure hemm_ASS_CHA_SUA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure hemm_ASS_CHA_SUA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure hemm_ASS_CHA_SUA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure hemm_ASS_CHA_SUA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure hemm_ASS_CHA_SUA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure hemm_ASS_CHA_SUA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure hemm_ASS_CHA_SUA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure hemm_ASS_CHA_SUA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure hemm_ASS_CHA_SUA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure hemm_ASS_CHA_SLA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure hemm_ASS_CHA_SLA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure hemm_ASS_CHA_SLA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure hemm_ASS_CHA_SLA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure hemm_ASS_CHA_SLA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure hemm_ASS_CHA_SLA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure hemm_ASS_CHA_SLA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure hemm_ASS_CHA_SLA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure hemm_ASS_CHA_SLA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure hemm_ASS_CHA_SLA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure hemm_ASS_CHA_SLA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure hemm_ASS_CHA_SLA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure hemm_ASS_CHA_SLA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure hemm_ASS_CHA_SLA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure hemm_ASS_CHA_SLA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CNB_ENABLED
#undef SFB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CNA_ENABLED 1
#define SFA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CHB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SUB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SUB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SUB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SUB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SUB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SUB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SUB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SUB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SUB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SUB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SUB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SUB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SUB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SUB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SUB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CHB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SLB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SLB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SLB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SLB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SLB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SLB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SLB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SLB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SLB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SLB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SLB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SLB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SLB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SLB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure hemm_ASS_CNA_SFA_CHB_SLB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CNA_ENABLED
#undef SFA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ASS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define EXP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CNB_ENABLED 1
#define SFB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure hemm_EXP_CHA_SUA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure hemm_EXP_CHA_SUA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure hemm_EXP_CHA_SUA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure hemm_EXP_CHA_SUA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure hemm_EXP_CHA_SUA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure hemm_EXP_CHA_SUA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure hemm_EXP_CHA_SUA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure hemm_EXP_CHA_SUA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure hemm_EXP_CHA_SUA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure hemm_EXP_CHA_SUA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure hemm_EXP_CHA_SUA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure hemm_EXP_CHA_SUA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure hemm_EXP_CHA_SUA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure hemm_EXP_CHA_SUA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure hemm_EXP_CHA_SUA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure hemm_EXP_CHA_SLA_CNB_SFB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure hemm_EXP_CHA_SLA_CNB_SFB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure hemm_EXP_CHA_SLA_CNB_SFB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure hemm_EXP_CHA_SLA_CNB_SFB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure hemm_EXP_CHA_SLA_CNB_SFB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure hemm_EXP_CHA_SLA_CNB_SFB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure hemm_EXP_CHA_SLA_CNB_SFB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure hemm_EXP_CHA_SLA_CNB_SFB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure hemm_EXP_CHA_SLA_CNB_SFB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure hemm_EXP_CHA_SLA_CNB_SFB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure hemm_EXP_CHA_SLA_CNB_SFB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure hemm_EXP_CHA_SLA_CNB_SFB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure hemm_EXP_CHA_SLA_CNB_SFB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure hemm_EXP_CHA_SLA_CNB_SFB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure hemm_EXP_CHA_SLA_CNB_SFB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CNB_ENABLED
#undef SFB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CNA_ENABLED 1
#define SFA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CHB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SUB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SUB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SUB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SUB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SUB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SUB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SUB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SUB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SUB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SUB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SUB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SUB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SUB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SUB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SUB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CHB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SLB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SLB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SLB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SLB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SLB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SLB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SLB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SLB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SLB_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SLB_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SLB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SLB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SLB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SLB_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure hemm_EXP_CNA_SFA_CHB_SLB_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CNA_ENABLED
#undef SFA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef EXP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef hemm_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define gemm_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SFA_ENABLED 1
#define SFB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ASS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TNB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TNA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TNB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TNB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TNB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TNB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TNB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TNB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TNB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TNB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TNB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TNB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TNB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TNB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TNB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TNB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TNB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TNA_ENABLED
#undef TNB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TSA_ENABLED 1
#define TNB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TNB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TNB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TNB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TNB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TNB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TNB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TNB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TNB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TNB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TNB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TNB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TNB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TNB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TNB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TNB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TSA_ENABLED
#undef TNB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define THA_ENABLED 1
#define TNB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TNB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TNB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TNB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TNB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TNB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TNB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TNB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TNB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TNB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TNB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TNB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TNB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TNB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TNB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TNB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef THA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TNB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TNA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TSB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TSB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TSB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TSB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TSB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TSB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TSB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TSB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TSB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TSB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TSB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TSB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TSB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TSB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_TSB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TNA_ENABLED
#undef TSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TSA_ENABLED 1
#define TSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TSB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TSB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TSB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TSB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TSB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TSB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TSB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TSB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TSB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TSB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TSB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TSB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TSB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TSB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_TSB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TSA_ENABLED
#undef TSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define THA_ENABLED 1
#define TSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TSB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TSB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TSB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TSB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TSB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TSB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TSB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TSB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TSB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TSB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TSB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TSB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TSB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TSB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_TSB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef THA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define THB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TNA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_THB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_THB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_THB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_THB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_THB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_THB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_THB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_THB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_THB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_THB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_THB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_THB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_THB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_THB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TNA_THB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TNA_ENABLED
#undef THB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TSA_ENABLED 1
#define THB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_THB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_THB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_THB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_THB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_THB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_THB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_THB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_THB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_THB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_THB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_THB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_THB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_THB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_THB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_TSA_THB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TSA_ENABLED
#undef THB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define THA_ENABLED 1
#define THB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_THB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_THB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_THB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_THB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_THB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_THB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_THB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_THB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_THB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_THB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_THB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_THB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_THB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_THB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_ASS_SFA_SFB_THA_THB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef THA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef THB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ASS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define EXP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TNB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TNA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TNB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TNB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TNB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TNB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TNB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TNB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TNB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TNB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TNB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TNB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TNB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TNB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TNB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TNB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TNB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TNA_ENABLED
#undef TNB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TSA_ENABLED 1
#define TNB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TNB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TNB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TNB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TNB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TNB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TNB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TNB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TNB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TNB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TNB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TNB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TNB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TNB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TNB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TNB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TSA_ENABLED
#undef TNB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define THA_ENABLED 1
#define TNB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TNB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TNB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TNB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TNB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TNB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TNB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TNB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TNB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TNB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TNB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TNB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TNB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TNB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TNB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TNB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef THA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TNB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TNA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TSB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TSB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TSB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TSB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TSB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TSB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TSB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TSB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TSB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TSB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TSB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TSB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TSB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TSB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_TSB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TNA_ENABLED
#undef TSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TSA_ENABLED 1
#define TSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TSB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TSB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TSB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TSB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TSB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TSB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TSB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TSB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TSB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TSB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TSB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TSB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TSB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TSB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_TSB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TSA_ENABLED
#undef TSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define THA_ENABLED 1
#define TSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TSB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TSB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TSB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TSB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TSB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TSB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TSB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TSB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TSB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TSB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TSB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TSB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TSB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TSB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_TSB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef THA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define THB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TNA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_THB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_THB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_THB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_THB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_THB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_THB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_THB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_THB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_THB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_THB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_THB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_THB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_THB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_THB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TNA_THB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TNA_ENABLED
#undef THB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TSA_ENABLED 1
#define THB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_THB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_THB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_THB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_THB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_THB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_THB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_THB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_THB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_THB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_THB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_THB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_THB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_THB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_THB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_TSA_THB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TSA_ENABLED
#undef THB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define THA_ENABLED 1
#define THB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_THB_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_THB_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_THB_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_THB_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_THB_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_THB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_THB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_THB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_THB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_THB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_THB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_THB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_THB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_THB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure gemm_EXP_SFA_SFB_THA_THB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulAdd@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef THA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef THB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef EXP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SFA_ENABLED
#undef SFB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef gemm_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatMulAdd_ENABLED

#undef CHECK_ASSERTION

end submodule routines