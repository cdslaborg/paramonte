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
!>  This file contains procedure implementations of [pm_matrixTrace](@ref pm_matrixTrace).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_matrixTrace) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_kind, only: RKD
    use pm_mathSqrt, only: getSqrt
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatTrace_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DEF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatTrace_DEF_XXX_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatTrace_DEF_XXX_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatTrace_DEF_XXX_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatTrace_DEF_XXX_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatTrace_DEF_XXX_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatTrace_DEF_XXX_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatTrace_DEF_XXX_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatTrace_DEF_XXX_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatTrace_DEF_XXX_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatTrace_DEF_XXX_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatTrace_DEF_XXX_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatTrace_DEF_XXX_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatTrace_DEF_XXX_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatTrace_DEF_XXX_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatTrace_DEF_XXX_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DEF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatTrace_RDP_XXX_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatTrace_RDP_XXX_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatTrace_RDP_XXX_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatTrace_RDP_XXX_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatTrace_RDP_XXX_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatTrace_RDP_XXX_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatTrace_RDP_XXX_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatTrace_RDP_XXX_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatTrace_RDP_XXX_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatTrace_RDP_XXX_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatTrace_RDP_XXX_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatTrace_RDP_XXX_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatTrace_RDP_XXX_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatTrace_RDP_XXX_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatTrace_RDP_XXX_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RFP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatTrace_RFP_UXD_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatTrace_RFP_UXD_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatTrace_RFP_UXD_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatTrace_RFP_UXD_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatTrace_RFP_UXD_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatTrace_RFP_UXD_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatTrace_RFP_UXD_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatTrace_RFP_UXD_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatTrace_RFP_UXD_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatTrace_RFP_UXD_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatTrace_RFP_UXD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatTrace_RFP_UXD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatTrace_RFP_UXD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatTrace_RFP_UXD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatTrace_RFP_UXD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatTrace_RFP_XLD_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatTrace_RFP_XLD_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatTrace_RFP_XLD_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatTrace_RFP_XLD_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatTrace_RFP_XLD_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatTrace_RFP_XLD_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatTrace_RFP_XLD_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatTrace_RFP_XLD_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatTrace_RFP_XLD_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatTrace_RFP_XLD_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatTrace_RFP_XLD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatTrace_RFP_XLD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatTrace_RFP_XLD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatTrace_RFP_XLD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatTrace_RFP_XLD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RFP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LFP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatTrace_LFP_UXD_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatTrace_LFP_UXD_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatTrace_LFP_UXD_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatTrace_LFP_UXD_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatTrace_LFP_UXD_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatTrace_LFP_UXD_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatTrace_LFP_UXD_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatTrace_LFP_UXD_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatTrace_LFP_UXD_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatTrace_LFP_UXD_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatTrace_LFP_UXD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatTrace_LFP_UXD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatTrace_LFP_UXD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatTrace_LFP_UXD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatTrace_LFP_UXD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatTrace_LFP_XLD_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatTrace_LFP_XLD_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatTrace_LFP_XLD_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatTrace_LFP_XLD_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatTrace_LFP_XLD_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatTrace_LFP_XLD_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatTrace_LFP_XLD_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatTrace_LFP_XLD_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatTrace_LFP_XLD_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatTrace_LFP_XLD_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatTrace_LFP_XLD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatTrace_LFP_XLD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatTrace_LFP_XLD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatTrace_LFP_XLD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatTrace_LFP_XLD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LFP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMatTrace_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatMulTrace_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DEF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatMulTrace_DEF_XXX_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatMulTrace_DEF_XXX_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatMulTrace_DEF_XXX_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatMulTrace_DEF_XXX_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatMulTrace_DEF_XXX_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatMulTrace_DEF_XXX_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatMulTrace_DEF_XXX_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatMulTrace_DEF_XXX_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatMulTrace_DEF_XXX_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatMulTrace_DEF_XXX_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatMulTrace_DEF_XXX_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatMulTrace_DEF_XXX_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatMulTrace_DEF_XXX_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatMulTrace_DEF_XXX_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatMulTrace_DEF_XXX_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DEF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatMulTrace_RDP_XXX_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatMulTrace_RDP_XXX_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatMulTrace_RDP_XXX_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatMulTrace_RDP_XXX_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatMulTrace_RDP_XXX_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatMulTrace_RDP_XXX_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatMulTrace_RDP_XXX_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatMulTrace_RDP_XXX_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatMulTrace_RDP_XXX_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatMulTrace_RDP_XXX_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatMulTrace_RDP_XXX_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatMulTrace_RDP_XXX_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatMulTrace_RDP_XXX_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatMulTrace_RDP_XXX_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatMulTrace_RDP_XXX_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RFP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatMulTrace_RFP_UXD_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatMulTrace_RFP_UXD_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatMulTrace_RFP_UXD_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatMulTrace_RFP_UXD_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatMulTrace_RFP_UXD_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatMulTrace_RFP_UXD_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatMulTrace_RFP_UXD_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatMulTrace_RFP_UXD_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatMulTrace_RFP_UXD_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatMulTrace_RFP_UXD_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatMulTrace_RFP_UXD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatMulTrace_RFP_UXD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatMulTrace_RFP_UXD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatMulTrace_RFP_UXD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatMulTrace_RFP_UXD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatMulTrace_RFP_XLD_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatMulTrace_RFP_XLD_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatMulTrace_RFP_XLD_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatMulTrace_RFP_XLD_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatMulTrace_RFP_XLD_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatMulTrace_RFP_XLD_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatMulTrace_RFP_XLD_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatMulTrace_RFP_XLD_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatMulTrace_RFP_XLD_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatMulTrace_RFP_XLD_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatMulTrace_RFP_XLD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatMulTrace_RFP_XLD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatMulTrace_RFP_XLD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatMulTrace_RFP_XLD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatMulTrace_RFP_XLD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RFP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LFP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatMulTrace_LFP_UXD_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatMulTrace_LFP_UXD_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatMulTrace_LFP_UXD_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatMulTrace_LFP_UXD_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatMulTrace_LFP_UXD_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatMulTrace_LFP_UXD_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatMulTrace_LFP_UXD_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatMulTrace_LFP_UXD_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatMulTrace_LFP_UXD_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatMulTrace_LFP_UXD_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatMulTrace_LFP_UXD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatMulTrace_LFP_UXD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatMulTrace_LFP_UXD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatMulTrace_LFP_UXD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatMulTrace_LFP_UXD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatMulTrace_LFP_XLD_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatMulTrace_LFP_XLD_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatMulTrace_LFP_XLD_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatMulTrace_LFP_XLD_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatMulTrace_LFP_XLD_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatMulTrace_LFP_XLD_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatMulTrace_LFP_XLD_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatMulTrace_LFP_XLD_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatMulTrace_LFP_XLD_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatMulTrace_LFP_XLD_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatMulTrace_LFP_XLD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatMulTrace_LFP_XLD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatMulTrace_LFP_XLD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatMulTrace_LFP_XLD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatMulTrace_LFP_XLD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LFP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMatMulTrace_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatMulTraceLog_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DEF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatMulTraceLog_DEF_XXX_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatMulTraceLog_DEF_XXX_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatMulTraceLog_DEF_XXX_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatMulTraceLog_DEF_XXX_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatMulTraceLog_DEF_XXX_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatMulTraceLog_DEF_XXX_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatMulTraceLog_DEF_XXX_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatMulTraceLog_DEF_XXX_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatMulTraceLog_DEF_XXX_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatMulTraceLog_DEF_XXX_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatMulTraceLog_DEF_XXX_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatMulTraceLog_DEF_XXX_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatMulTraceLog_DEF_XXX_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatMulTraceLog_DEF_XXX_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatMulTraceLog_DEF_XXX_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DEF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatMulTraceLog_RDP_XXX_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatMulTraceLog_RDP_XXX_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatMulTraceLog_RDP_XXX_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatMulTraceLog_RDP_XXX_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatMulTraceLog_RDP_XXX_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatMulTraceLog_RDP_XXX_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatMulTraceLog_RDP_XXX_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatMulTraceLog_RDP_XXX_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatMulTraceLog_RDP_XXX_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatMulTraceLog_RDP_XXX_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatMulTraceLog_RDP_XXX_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatMulTraceLog_RDP_XXX_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatMulTraceLog_RDP_XXX_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatMulTraceLog_RDP_XXX_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatMulTraceLog_RDP_XXX_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RFP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatMulTraceLog_RFP_UXD_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatMulTraceLog_RFP_UXD_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatMulTraceLog_RFP_UXD_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatMulTraceLog_RFP_UXD_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatMulTraceLog_RFP_UXD_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatMulTraceLog_RFP_UXD_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatMulTraceLog_RFP_UXD_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatMulTraceLog_RFP_UXD_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatMulTraceLog_RFP_UXD_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatMulTraceLog_RFP_UXD_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatMulTraceLog_RFP_UXD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatMulTraceLog_RFP_UXD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatMulTraceLog_RFP_UXD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatMulTraceLog_RFP_UXD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatMulTraceLog_RFP_UXD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatMulTraceLog_RFP_XLD_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatMulTraceLog_RFP_XLD_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatMulTraceLog_RFP_XLD_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatMulTraceLog_RFP_XLD_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatMulTraceLog_RFP_XLD_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatMulTraceLog_RFP_XLD_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatMulTraceLog_RFP_XLD_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatMulTraceLog_RFP_XLD_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatMulTraceLog_RFP_XLD_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatMulTraceLog_RFP_XLD_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatMulTraceLog_RFP_XLD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatMulTraceLog_RFP_XLD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatMulTraceLog_RFP_XLD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatMulTraceLog_RFP_XLD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatMulTraceLog_RFP_XLD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RFP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LFP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatMulTraceLog_LFP_UXD_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatMulTraceLog_LFP_UXD_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatMulTraceLog_LFP_UXD_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatMulTraceLog_LFP_UXD_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatMulTraceLog_LFP_UXD_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatMulTraceLog_LFP_UXD_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatMulTraceLog_LFP_UXD_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatMulTraceLog_LFP_UXD_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatMulTraceLog_LFP_UXD_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatMulTraceLog_LFP_UXD_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatMulTraceLog_LFP_UXD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatMulTraceLog_LFP_UXD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatMulTraceLog_LFP_UXD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatMulTraceLog_LFP_UXD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatMulTraceLog_LFP_UXD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMatMulTraceLog_LFP_XLD_IK5
        use pm_kind, only: IKG => IK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMatMulTraceLog_LFP_XLD_IK4
        use pm_kind, only: IKG => IK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMatMulTraceLog_LFP_XLD_IK3
        use pm_kind, only: IKG => IK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMatMulTraceLog_LFP_XLD_IK2
        use pm_kind, only: IKG => IK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMatMulTraceLog_LFP_XLD_IK1
        use pm_kind, only: IKG => IK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatMulTraceLog_LFP_XLD_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatMulTraceLog_LFP_XLD_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatMulTraceLog_LFP_XLD_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatMulTraceLog_LFP_XLD_CK2
        use pm_kind, only: CKG => CK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatMulTraceLog_LFP_XLD_CK1
        use pm_kind, only: CKG => CK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatMulTraceLog_LFP_XLD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatMulTraceLog_LFP_XLD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatMulTraceLog_LFP_XLD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatMulTraceLog_LFP_XLD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatMulTraceLog_LFP_XLD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_matrixTrace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LFP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMatMulTraceLog_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines