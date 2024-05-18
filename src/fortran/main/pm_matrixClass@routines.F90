
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
!>  This file contains procedure implementations of [pm_matrixClass](@ref pm_matrixClass).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_matrixClass) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_matrixPack, only: isMatPack
    use pm_matrixCopy, only: setMatCopy
    use pm_matrixChol, only: setMatChol
    use pm_matrixChol, only: recursion
    use pm_matrixChol, only: nothing
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isMatClass_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PosDef_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ful_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isMatClassPosDefFulRDP_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isMatClassPosDefFulRDP_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isMatClassPosDefFulRDP_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isMatClassPosDefFulRDP_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isMatClassPosDefFulRDP_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isMatClassPosDefFulRDP_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isMatClassPosDefFulRDP_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isMatClassPosDefFulRDP_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isMatClassPosDefFulRDP_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isMatClassPosDefFulRDP_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ful_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Upp_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isMatClassPosDefUppRDP_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isMatClassPosDefUppRDP_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isMatClassPosDefUppRDP_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isMatClassPosDefUppRDP_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isMatClassPosDefUppRDP_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isMatClassPosDefUppRDP_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isMatClassPosDefUppRDP_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isMatClassPosDefUppRDP_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isMatClassPosDefUppRDP_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isMatClassPosDefUppRDP_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Upp_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Low_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isMatClassPosDefLowRDP_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isMatClassPosDefLowRDP_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isMatClassPosDefLowRDP_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isMatClassPosDefLowRDP_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isMatClassPosDefLowRDP_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isMatClassPosDefLowRDP_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isMatClassPosDefLowRDP_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isMatClassPosDefLowRDP_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isMatClassPosDefLowRDP_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isMatClassPosDefLowRDP_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Low_ENABLED

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

#define Upp_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isMatClassPosDefUppRFP_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isMatClassPosDefUppRFP_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isMatClassPosDefUppRFP_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isMatClassPosDefUppRFP_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isMatClassPosDefUppRFP_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isMatClassPosDefUppRFP_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isMatClassPosDefUppRFP_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isMatClassPosDefUppRFP_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isMatClassPosDefUppRFP_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isMatClassPosDefUppRFP_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Upp_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Low_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isMatClassPosDefLowRFP_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isMatClassPosDefLowRFP_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isMatClassPosDefLowRFP_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isMatClassPosDefLowRFP_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isMatClassPosDefLowRFP_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isMatClassPosDefLowRFP_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isMatClassPosDefLowRFP_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isMatClassPosDefLowRFP_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isMatClassPosDefLowRFP_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isMatClassPosDefLowRFP_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Low_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RFP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PosDef_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Symm_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isMatClassSymm_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isMatClassSymm_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isMatClassSymm_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isMatClassSymm_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isMatClassSymm_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure isMatClassSymm_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure isMatClassSymm_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure isMatClassSymm_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure isMatClassSymm_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure isMatClassSymm_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure isMatClassSymm_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure isMatClassSymm_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure isMatClassSymm_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure isMatClassSymm_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure isMatClassSymm_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isMatClassSymm_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isMatClassSymm_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isMatClassSymm_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isMatClassSymm_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isMatClassSymm_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isMatClassSymm_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isMatClassSymm_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isMatClassSymm_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isMatClassSymm_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isMatClassSymm_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Symm_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Herm_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isMatClassHerm_SK5
        use pm_kind, only: SKC => SK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isMatClassHerm_SK4
        use pm_kind, only: SKC => SK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isMatClassHerm_SK3
        use pm_kind, only: SKC => SK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isMatClassHerm_SK2
        use pm_kind, only: SKC => SK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isMatClassHerm_SK1
        use pm_kind, only: SKC => SK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure isMatClassHerm_IK5
        use pm_kind, only: IKC => IK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure isMatClassHerm_IK4
        use pm_kind, only: IKC => IK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure isMatClassHerm_IK3
        use pm_kind, only: IKC => IK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure isMatClassHerm_IK2
        use pm_kind, only: IKC => IK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure isMatClassHerm_IK1
        use pm_kind, only: IKC => IK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure isMatClassHerm_LK5
        use pm_kind, only: LKC => LK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure isMatClassHerm_LK4
        use pm_kind, only: LKC => LK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure isMatClassHerm_LK3
        use pm_kind, only: LKC => LK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure isMatClassHerm_LK2
        use pm_kind, only: LKC => LK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure isMatClassHerm_LK1
        use pm_kind, only: LKC => LK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isMatClassHerm_CK5
        use pm_kind, only: CKC => CK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isMatClassHerm_CK4
        use pm_kind, only: CKC => CK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isMatClassHerm_CK3
        use pm_kind, only: CKC => CK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isMatClassHerm_CK2
        use pm_kind, only: CKC => CK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isMatClassHerm_CK1
        use pm_kind, only: CKC => CK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isMatClassHerm_RK5
        use pm_kind, only: RKC => RK5
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isMatClassHerm_RK4
        use pm_kind, only: RKC => RK4
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isMatClassHerm_RK3
        use pm_kind, only: RKC => RK3
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isMatClassHerm_RK2
        use pm_kind, only: RKC => RK2
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isMatClassHerm_RK1
        use pm_kind, only: RKC => RK1
#include "pm_matrixClass@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Herm_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isMatClass_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
