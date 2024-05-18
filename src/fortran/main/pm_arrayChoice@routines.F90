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
!>  This file contains procedure implementations of [pm_arrayChoice](@ref pm_arrayChoice).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_arrayChoice) routines ! LCOV_EXCL_LINE

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
    use pm_distUnif, only: rngf
    use pm_arrayRange, only: setRange
    use pm_arraySearch, only: getBin
    use pm_arrayShuffle, only: setShuffled
    use pm_arraySpace, only: setLinSpace
    use pm_distUnif, only: setUnifRand
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getChoice_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getChoiceRNGD_D0_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getChoiceRNGD_D0_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getChoiceRNGD_D0_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getChoiceRNGD_D0_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getChoiceRNGD_D0_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_S1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getChoiceRNGD_D0_S1_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getChoiceRNGD_D0_S1_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getChoiceRNGD_D0_S1_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getChoiceRNGD_D0_S1_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getChoiceRNGD_D0_S1_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_S1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getChoiceRNGD_D1_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getChoiceRNGD_D1_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getChoiceRNGD_D1_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getChoiceRNGD_D1_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getChoiceRNGD_D1_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getChoiceRNGD_D1_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getChoiceRNGD_D1_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getChoiceRNGD_D1_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getChoiceRNGD_D1_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getChoiceRNGD_D1_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getChoiceRNGD_D1_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getChoiceRNGD_D1_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getChoiceRNGD_D1_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getChoiceRNGD_D1_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getChoiceRNGD_D1_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getChoiceRNGD_D1_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getChoiceRNGD_D1_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getChoiceRNGD_D1_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getChoiceRNGD_D1_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getChoiceRNGD_D1_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getChoiceRNGD_D1_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getChoiceRNGD_D1_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getChoiceRNGD_D1_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getChoiceRNGD_D1_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getChoiceRNGD_D1_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getChoiceRNGD_D1_D1_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getChoiceRNGD_D1_D1_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getChoiceRNGD_D1_D1_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getChoiceRNGD_D1_D1_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getChoiceRNGD_D1_D1_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getChoiceRNGD_D1_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getChoiceRNGD_D1_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getChoiceRNGD_D1_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getChoiceRNGD_D1_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getChoiceRNGD_D1_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getChoiceRNGD_D1_D1_LK5
        use pm_kind, only: LKC => LK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getChoiceRNGD_D1_D1_LK4
        use pm_kind, only: LKC => LK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getChoiceRNGD_D1_D1_LK3
        use pm_kind, only: LKC => LK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getChoiceRNGD_D1_D1_LK2
        use pm_kind, only: LKC => LK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getChoiceRNGD_D1_D1_LK1
        use pm_kind, only: LKC => LK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getChoiceRNGD_D1_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getChoiceRNGD_D1_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getChoiceRNGD_D1_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getChoiceRNGD_D1_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getChoiceRNGD_D1_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getChoiceRNGD_D1_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getChoiceRNGD_D1_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getChoiceRNGD_D1_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getChoiceRNGD_D1_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getChoiceRNGD_D1_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getChoice_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setChoice_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setChoiceRNGF_D0_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setChoiceRNGF_D0_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setChoiceRNGF_D0_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setChoiceRNGF_D0_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setChoiceRNGF_D0_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arrayChoice@routines.inc.F90"
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
    module procedure setChoiceRNGF_D1_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setChoiceRNGF_D1_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setChoiceRNGF_D1_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setChoiceRNGF_D1_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setChoiceRNGF_D1_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setChoiceRNGF_D1_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setChoiceRNGF_D1_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setChoiceRNGF_D1_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setChoiceRNGF_D1_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setChoiceRNGF_D1_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setChoiceRNGF_D1_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setChoiceRNGF_D1_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setChoiceRNGF_D1_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setChoiceRNGF_D1_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setChoiceRNGF_D1_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setChoiceRNGF_D1_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setChoiceRNGF_D1_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setChoiceRNGF_D1_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setChoiceRNGF_D1_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setChoiceRNGF_D1_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setChoiceRNGF_D1_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setChoiceRNGF_D1_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setChoiceRNGF_D1_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setChoiceRNGF_D1_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setChoiceRNGF_D1_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setChoiceRNGF_D1_D1_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setChoiceRNGF_D1_D1_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setChoiceRNGF_D1_D1_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setChoiceRNGF_D1_D1_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setChoiceRNGF_D1_D1_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setChoiceRNGF_D1_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setChoiceRNGF_D1_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setChoiceRNGF_D1_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setChoiceRNGF_D1_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setChoiceRNGF_D1_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setChoiceRNGF_D1_D1_LK5
        use pm_kind, only: LKC => LK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setChoiceRNGF_D1_D1_LK4
        use pm_kind, only: LKC => LK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setChoiceRNGF_D1_D1_LK3
        use pm_kind, only: LKC => LK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setChoiceRNGF_D1_D1_LK2
        use pm_kind, only: LKC => LK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setChoiceRNGF_D1_D1_LK1
        use pm_kind, only: LKC => LK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setChoiceRNGF_D1_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setChoiceRNGF_D1_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setChoiceRNGF_D1_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setChoiceRNGF_D1_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setChoiceRNGF_D1_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setChoiceRNGF_D1_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setChoiceRNGF_D1_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setChoiceRNGF_D1_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setChoiceRNGF_D1_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setChoiceRNGF_D1_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setChoiceRNGX_D0_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setChoiceRNGX_D0_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setChoiceRNGX_D0_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setChoiceRNGX_D0_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setChoiceRNGX_D0_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arrayChoice@routines.inc.F90"
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
    module procedure setChoiceRNGX_D1_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setChoiceRNGX_D1_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setChoiceRNGX_D1_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setChoiceRNGX_D1_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setChoiceRNGX_D1_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setChoiceRNGX_D1_D0_IK5
        use pm_kind, only: IKC => IK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setChoiceRNGX_D1_D0_IK4
        use pm_kind, only: IKC => IK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setChoiceRNGX_D1_D0_IK3
        use pm_kind, only: IKC => IK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setChoiceRNGX_D1_D0_IK2
        use pm_kind, only: IKC => IK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setChoiceRNGX_D1_D0_IK1
        use pm_kind, only: IKC => IK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setChoiceRNGX_D1_D0_LK5
        use pm_kind, only: LKC => LK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setChoiceRNGX_D1_D0_LK4
        use pm_kind, only: LKC => LK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setChoiceRNGX_D1_D0_LK3
        use pm_kind, only: LKC => LK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setChoiceRNGX_D1_D0_LK2
        use pm_kind, only: LKC => LK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setChoiceRNGX_D1_D0_LK1
        use pm_kind, only: LKC => LK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setChoiceRNGX_D1_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setChoiceRNGX_D1_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setChoiceRNGX_D1_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setChoiceRNGX_D1_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setChoiceRNGX_D1_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setChoiceRNGX_D1_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setChoiceRNGX_D1_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setChoiceRNGX_D1_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setChoiceRNGX_D1_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setChoiceRNGX_D1_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setChoiceRNGX_D1_D1_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setChoiceRNGX_D1_D1_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setChoiceRNGX_D1_D1_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setChoiceRNGX_D1_D1_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setChoiceRNGX_D1_D1_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setChoiceRNGX_D1_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setChoiceRNGX_D1_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setChoiceRNGX_D1_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setChoiceRNGX_D1_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setChoiceRNGX_D1_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setChoiceRNGX_D1_D1_LK5
        use pm_kind, only: LKC => LK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setChoiceRNGX_D1_D1_LK4
        use pm_kind, only: LKC => LK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setChoiceRNGX_D1_D1_LK3
        use pm_kind, only: LKC => LK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setChoiceRNGX_D1_D1_LK2
        use pm_kind, only: LKC => LK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setChoiceRNGX_D1_D1_LK1
        use pm_kind, only: LKC => LK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setChoiceRNGX_D1_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setChoiceRNGX_D1_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setChoiceRNGX_D1_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setChoiceRNGX_D1_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setChoiceRNGX_D1_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setChoiceRNGX_D1_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setChoiceRNGX_D1_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setChoiceRNGX_D1_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setChoiceRNGX_D1_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setChoiceRNGX_D1_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_arrayChoice@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setChoice_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines
