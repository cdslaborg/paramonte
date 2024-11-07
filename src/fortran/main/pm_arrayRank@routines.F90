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
!>  This file contains procedure implementations of [pm_arrayRankFrac](@ref pm_arrayRankFrac).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_arrayRank) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_arraySort, only: setSorted
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Dense_ENABLED 1

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
    module procedure getRankDenseDefCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankDenseDefCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankDenseDefCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankDenseDefCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankDenseDefCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure getRankDenseDefCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankDenseDefCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankDenseDefCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankDenseDefCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankDenseDefCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRankDenseDefCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRankDenseDefCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRankDenseDefCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRankDenseDefCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRankDenseDefCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getRankDenseDefCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getRankDenseDefCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getRankDenseDefCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getRankDenseDefCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getRankDenseDefCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getRankDenseDefCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getRankDenseDefCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getRankDenseDefCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getRankDenseDefCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getRankDenseDefCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRankDenseDefCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRankDenseDefCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRankDenseDefCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRankDenseDefCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRankDenseDefCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getRankDenseDefCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankDenseDefCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankDenseDefCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankDenseDefCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankDenseDefCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getRankDenseDefCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure getRankDenseCusCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankDenseCusCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankDenseCusCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankDenseCusCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankDenseCusCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRankDenseCusCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankDenseCusCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankDenseCusCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankDenseCusCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankDenseCusCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRankDenseCusCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRankDenseCusCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRankDenseCusCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRankDenseCusCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRankDenseCusCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getRankDenseCusCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getRankDenseCusCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getRankDenseCusCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getRankDenseCusCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getRankDenseCusCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getRankDenseCusCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getRankDenseCusCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getRankDenseCusCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getRankDenseCusCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getRankDenseCusCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRankDenseCusCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRankDenseCusCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRankDenseCusCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRankDenseCusCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRankDenseCusCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getRankDenseCusCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankDenseCusCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankDenseCusCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankDenseCusCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankDenseCusCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getRankDenseCusCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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

#undef Dense_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Dense_ENABLED 1

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
    module procedure setRankDenseDefCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankDenseDefCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankDenseDefCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankDenseDefCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankDenseDefCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure setRankDenseDefCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankDenseDefCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankDenseDefCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankDenseDefCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankDenseDefCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRankDenseDefCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRankDenseDefCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRankDenseDefCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRankDenseDefCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRankDenseDefCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setRankDenseDefCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setRankDenseDefCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setRankDenseDefCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setRankDenseDefCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setRankDenseDefCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setRankDenseDefCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setRankDenseDefCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setRankDenseDefCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setRankDenseDefCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setRankDenseDefCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRankDenseDefCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRankDenseDefCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRankDenseDefCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRankDenseDefCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRankDenseDefCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRankDenseDefCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankDenseDefCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankDenseDefCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankDenseDefCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankDenseDefCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRankDenseDefCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure setRankDenseCusCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankDenseCusCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankDenseCusCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankDenseCusCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankDenseCusCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRankDenseCusCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankDenseCusCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankDenseCusCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankDenseCusCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankDenseCusCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRankDenseCusCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRankDenseCusCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRankDenseCusCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRankDenseCusCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRankDenseCusCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setRankDenseCusCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setRankDenseCusCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setRankDenseCusCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setRankDenseCusCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setRankDenseCusCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setRankDenseCusCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setRankDenseCusCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setRankDenseCusCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setRankDenseCusCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setRankDenseCusCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRankDenseCusCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRankDenseCusCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRankDenseCusCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRankDenseCusCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRankDenseCusCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRankDenseCusCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankDenseCusCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankDenseCusCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankDenseCusCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankDenseCusCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRankDenseCusCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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

#undef Dense_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Fractional_ENABLED 1

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
    module procedure getRankFractionalDefCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankFractionalDefCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankFractionalDefCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankFractionalDefCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankFractionalDefCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure getRankFractionalDefCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankFractionalDefCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankFractionalDefCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankFractionalDefCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankFractionalDefCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRankFractionalDefCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRankFractionalDefCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRankFractionalDefCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRankFractionalDefCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRankFractionalDefCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getRankFractionalDefCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getRankFractionalDefCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getRankFractionalDefCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getRankFractionalDefCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getRankFractionalDefCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getRankFractionalDefCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getRankFractionalDefCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getRankFractionalDefCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getRankFractionalDefCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getRankFractionalDefCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRankFractionalDefCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRankFractionalDefCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRankFractionalDefCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRankFractionalDefCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRankFractionalDefCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getRankFractionalDefCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankFractionalDefCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankFractionalDefCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankFractionalDefCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankFractionalDefCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getRankFractionalDefCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure getRankFractionalCusCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankFractionalCusCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankFractionalCusCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankFractionalCusCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankFractionalCusCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRankFractionalCusCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankFractionalCusCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankFractionalCusCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankFractionalCusCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankFractionalCusCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRankFractionalCusCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRankFractionalCusCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRankFractionalCusCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRankFractionalCusCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRankFractionalCusCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getRankFractionalCusCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getRankFractionalCusCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getRankFractionalCusCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getRankFractionalCusCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getRankFractionalCusCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getRankFractionalCusCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getRankFractionalCusCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getRankFractionalCusCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getRankFractionalCusCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getRankFractionalCusCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRankFractionalCusCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRankFractionalCusCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRankFractionalCusCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRankFractionalCusCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRankFractionalCusCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getRankFractionalCusCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankFractionalCusCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankFractionalCusCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankFractionalCusCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankFractionalCusCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getRankFractionalCusCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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

#undef Fractional_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Fractional_ENABLED 1

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
    module procedure setRankFractionalDefCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankFractionalDefCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankFractionalDefCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankFractionalDefCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankFractionalDefCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure setRankFractionalDefCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankFractionalDefCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankFractionalDefCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankFractionalDefCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankFractionalDefCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRankFractionalDefCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRankFractionalDefCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRankFractionalDefCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRankFractionalDefCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRankFractionalDefCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setRankFractionalDefCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setRankFractionalDefCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setRankFractionalDefCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setRankFractionalDefCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setRankFractionalDefCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setRankFractionalDefCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setRankFractionalDefCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setRankFractionalDefCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setRankFractionalDefCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setRankFractionalDefCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRankFractionalDefCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRankFractionalDefCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRankFractionalDefCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRankFractionalDefCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRankFractionalDefCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRankFractionalDefCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankFractionalDefCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankFractionalDefCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankFractionalDefCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankFractionalDefCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRankFractionalDefCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure setRankFractionalCusCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankFractionalCusCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankFractionalCusCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankFractionalCusCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankFractionalCusCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRankFractionalCusCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankFractionalCusCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankFractionalCusCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankFractionalCusCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankFractionalCusCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRankFractionalCusCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRankFractionalCusCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRankFractionalCusCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRankFractionalCusCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRankFractionalCusCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setRankFractionalCusCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setRankFractionalCusCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setRankFractionalCusCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setRankFractionalCusCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setRankFractionalCusCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setRankFractionalCusCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setRankFractionalCusCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setRankFractionalCusCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setRankFractionalCusCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setRankFractionalCusCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRankFractionalCusCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRankFractionalCusCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRankFractionalCusCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRankFractionalCusCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRankFractionalCusCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRankFractionalCusCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankFractionalCusCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankFractionalCusCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankFractionalCusCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankFractionalCusCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRankFractionalCusCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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

#undef Fractional_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Modified_ENABLED 1

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
    module procedure getRankModifiedDefCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankModifiedDefCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankModifiedDefCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankModifiedDefCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankModifiedDefCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure getRankModifiedDefCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankModifiedDefCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankModifiedDefCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankModifiedDefCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankModifiedDefCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRankModifiedDefCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRankModifiedDefCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRankModifiedDefCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRankModifiedDefCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRankModifiedDefCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getRankModifiedDefCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getRankModifiedDefCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getRankModifiedDefCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getRankModifiedDefCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getRankModifiedDefCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getRankModifiedDefCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getRankModifiedDefCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getRankModifiedDefCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getRankModifiedDefCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getRankModifiedDefCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRankModifiedDefCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRankModifiedDefCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRankModifiedDefCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRankModifiedDefCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRankModifiedDefCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getRankModifiedDefCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankModifiedDefCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankModifiedDefCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankModifiedDefCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankModifiedDefCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getRankModifiedDefCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure getRankModifiedCusCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankModifiedCusCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankModifiedCusCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankModifiedCusCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankModifiedCusCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRankModifiedCusCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankModifiedCusCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankModifiedCusCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankModifiedCusCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankModifiedCusCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRankModifiedCusCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRankModifiedCusCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRankModifiedCusCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRankModifiedCusCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRankModifiedCusCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getRankModifiedCusCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getRankModifiedCusCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getRankModifiedCusCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getRankModifiedCusCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getRankModifiedCusCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getRankModifiedCusCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getRankModifiedCusCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getRankModifiedCusCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getRankModifiedCusCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getRankModifiedCusCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRankModifiedCusCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRankModifiedCusCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRankModifiedCusCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRankModifiedCusCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRankModifiedCusCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getRankModifiedCusCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankModifiedCusCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankModifiedCusCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankModifiedCusCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankModifiedCusCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getRankModifiedCusCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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

#undef Modified_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Modified_ENABLED 1

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
    module procedure setRankModifiedDefCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankModifiedDefCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankModifiedDefCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankModifiedDefCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankModifiedDefCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure setRankModifiedDefCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankModifiedDefCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankModifiedDefCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankModifiedDefCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankModifiedDefCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRankModifiedDefCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRankModifiedDefCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRankModifiedDefCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRankModifiedDefCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRankModifiedDefCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setRankModifiedDefCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setRankModifiedDefCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setRankModifiedDefCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setRankModifiedDefCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setRankModifiedDefCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setRankModifiedDefCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setRankModifiedDefCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setRankModifiedDefCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setRankModifiedDefCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setRankModifiedDefCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRankModifiedDefCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRankModifiedDefCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRankModifiedDefCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRankModifiedDefCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRankModifiedDefCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRankModifiedDefCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankModifiedDefCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankModifiedDefCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankModifiedDefCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankModifiedDefCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRankModifiedDefCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure setRankModifiedCusCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankModifiedCusCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankModifiedCusCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankModifiedCusCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankModifiedCusCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRankModifiedCusCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankModifiedCusCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankModifiedCusCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankModifiedCusCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankModifiedCusCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRankModifiedCusCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRankModifiedCusCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRankModifiedCusCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRankModifiedCusCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRankModifiedCusCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setRankModifiedCusCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setRankModifiedCusCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setRankModifiedCusCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setRankModifiedCusCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setRankModifiedCusCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setRankModifiedCusCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setRankModifiedCusCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setRankModifiedCusCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setRankModifiedCusCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setRankModifiedCusCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRankModifiedCusCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRankModifiedCusCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRankModifiedCusCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRankModifiedCusCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRankModifiedCusCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRankModifiedCusCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankModifiedCusCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankModifiedCusCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankModifiedCusCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankModifiedCusCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRankModifiedCusCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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

#undef Modified_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ordinal_ENABLED 1

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
    module procedure getRankOrdinalDefCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankOrdinalDefCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankOrdinalDefCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankOrdinalDefCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankOrdinalDefCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure getRankOrdinalDefCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankOrdinalDefCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankOrdinalDefCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankOrdinalDefCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankOrdinalDefCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRankOrdinalDefCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRankOrdinalDefCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRankOrdinalDefCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRankOrdinalDefCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRankOrdinalDefCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getRankOrdinalDefCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getRankOrdinalDefCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getRankOrdinalDefCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getRankOrdinalDefCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getRankOrdinalDefCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getRankOrdinalDefCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getRankOrdinalDefCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getRankOrdinalDefCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getRankOrdinalDefCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getRankOrdinalDefCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRankOrdinalDefCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRankOrdinalDefCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRankOrdinalDefCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRankOrdinalDefCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRankOrdinalDefCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getRankOrdinalDefCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankOrdinalDefCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankOrdinalDefCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankOrdinalDefCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankOrdinalDefCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getRankOrdinalDefCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure getRankOrdinalCusCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankOrdinalCusCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankOrdinalCusCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankOrdinalCusCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankOrdinalCusCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRankOrdinalCusCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankOrdinalCusCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankOrdinalCusCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankOrdinalCusCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankOrdinalCusCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRankOrdinalCusCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRankOrdinalCusCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRankOrdinalCusCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRankOrdinalCusCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRankOrdinalCusCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getRankOrdinalCusCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getRankOrdinalCusCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getRankOrdinalCusCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getRankOrdinalCusCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getRankOrdinalCusCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getRankOrdinalCusCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getRankOrdinalCusCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getRankOrdinalCusCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getRankOrdinalCusCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getRankOrdinalCusCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRankOrdinalCusCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRankOrdinalCusCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRankOrdinalCusCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRankOrdinalCusCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRankOrdinalCusCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getRankOrdinalCusCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankOrdinalCusCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankOrdinalCusCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankOrdinalCusCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankOrdinalCusCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getRankOrdinalCusCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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

#undef Ordinal_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ordinal_ENABLED 1

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
    module procedure setRankOrdinalDefCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankOrdinalDefCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankOrdinalDefCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankOrdinalDefCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankOrdinalDefCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure setRankOrdinalDefCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankOrdinalDefCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankOrdinalDefCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankOrdinalDefCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankOrdinalDefCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRankOrdinalDefCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRankOrdinalDefCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRankOrdinalDefCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRankOrdinalDefCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRankOrdinalDefCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setRankOrdinalDefCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setRankOrdinalDefCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setRankOrdinalDefCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setRankOrdinalDefCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setRankOrdinalDefCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setRankOrdinalDefCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setRankOrdinalDefCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setRankOrdinalDefCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setRankOrdinalDefCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setRankOrdinalDefCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRankOrdinalDefCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRankOrdinalDefCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRankOrdinalDefCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRankOrdinalDefCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRankOrdinalDefCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRankOrdinalDefCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankOrdinalDefCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankOrdinalDefCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankOrdinalDefCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankOrdinalDefCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRankOrdinalDefCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure setRankOrdinalCusCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankOrdinalCusCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankOrdinalCusCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankOrdinalCusCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankOrdinalCusCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRankOrdinalCusCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankOrdinalCusCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankOrdinalCusCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankOrdinalCusCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankOrdinalCusCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRankOrdinalCusCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRankOrdinalCusCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRankOrdinalCusCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRankOrdinalCusCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRankOrdinalCusCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setRankOrdinalCusCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setRankOrdinalCusCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setRankOrdinalCusCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setRankOrdinalCusCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setRankOrdinalCusCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setRankOrdinalCusCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setRankOrdinalCusCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setRankOrdinalCusCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setRankOrdinalCusCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setRankOrdinalCusCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRankOrdinalCusCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRankOrdinalCusCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRankOrdinalCusCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRankOrdinalCusCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRankOrdinalCusCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRankOrdinalCusCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankOrdinalCusCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankOrdinalCusCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankOrdinalCusCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankOrdinalCusCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRankOrdinalCusCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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

#undef Ordinal_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Standard_ENABLED 1

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
    module procedure getRankStandardDefCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankStandardDefCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankStandardDefCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankStandardDefCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankStandardDefCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure getRankStandardDefCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankStandardDefCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankStandardDefCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankStandardDefCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankStandardDefCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRankStandardDefCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRankStandardDefCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRankStandardDefCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRankStandardDefCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRankStandardDefCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getRankStandardDefCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getRankStandardDefCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getRankStandardDefCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getRankStandardDefCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getRankStandardDefCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getRankStandardDefCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getRankStandardDefCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getRankStandardDefCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getRankStandardDefCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getRankStandardDefCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRankStandardDefCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRankStandardDefCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRankStandardDefCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRankStandardDefCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRankStandardDefCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getRankStandardDefCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankStandardDefCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankStandardDefCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankStandardDefCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankStandardDefCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getRankStandardDefCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure getRankStandardCusCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankStandardCusCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankStandardCusCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankStandardCusCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankStandardCusCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRankStandardCusCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankStandardCusCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankStandardCusCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankStandardCusCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankStandardCusCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRankStandardCusCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRankStandardCusCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRankStandardCusCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRankStandardCusCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRankStandardCusCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getRankStandardCusCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getRankStandardCusCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getRankStandardCusCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getRankStandardCusCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getRankStandardCusCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getRankStandardCusCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getRankStandardCusCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getRankStandardCusCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getRankStandardCusCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getRankStandardCusCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRankStandardCusCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRankStandardCusCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRankStandardCusCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRankStandardCusCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRankStandardCusCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getRankStandardCusCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRankStandardCusCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRankStandardCusCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRankStandardCusCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRankStandardCusCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getRankStandardCusCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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

#undef Standard_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Standard_ENABLED 1

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
    module procedure setRankStandardDefCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankStandardDefCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankStandardDefCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankStandardDefCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankStandardDefCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure setRankStandardDefCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankStandardDefCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankStandardDefCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankStandardDefCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankStandardDefCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRankStandardDefCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRankStandardDefCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRankStandardDefCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRankStandardDefCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRankStandardDefCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setRankStandardDefCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setRankStandardDefCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setRankStandardDefCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setRankStandardDefCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setRankStandardDefCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setRankStandardDefCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setRankStandardDefCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setRankStandardDefCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setRankStandardDefCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setRankStandardDefCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRankStandardDefCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRankStandardDefCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRankStandardDefCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRankStandardDefCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRankStandardDefCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRankStandardDefCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankStandardDefCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankStandardDefCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankStandardDefCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankStandardDefCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRankStandardDefCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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
    module procedure setRankStandardCusCom_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankStandardCusCom_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankStandardCusCom_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankStandardCusCom_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankStandardCusCom_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRankStandardCusCom_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankStandardCusCom_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankStandardCusCom_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankStandardCusCom_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankStandardCusCom_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRankStandardCusCom_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRankStandardCusCom_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRankStandardCusCom_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRankStandardCusCom_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRankStandardCusCom_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setRankStandardCusCom_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setRankStandardCusCom_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setRankStandardCusCom_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setRankStandardCusCom_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setRankStandardCusCom_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setRankStandardCusCom_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setRankStandardCusCom_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setRankStandardCusCom_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setRankStandardCusCom_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setRankStandardCusCom_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRankStandardCusCom_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRankStandardCusCom_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRankStandardCusCom_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRankStandardCusCom_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRankStandardCusCom_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRankStandardCusCom_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRankStandardCusCom_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRankStandardCusCom_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRankStandardCusCom_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRankStandardCusCom_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRankStandardCusCom_D1_BSSK
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
#include "pm_arrayRank@routines.inc.F90"
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

#undef Standard_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines