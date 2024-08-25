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
!>  This file contains procedure implementations of [pm_mathSum](@ref pm_mathSum).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, August 8, 2024, 10:23 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_mathSum) routines ! LCOV_EXCL_LINE

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

#define getSum_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getSumDef_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getSumDef_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getSumDef_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getSumDef_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getSumDef_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getSumDef_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getSumDef_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getSumDef_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getSumDef_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getSumDef_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ite_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getSumIte_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getSumIte_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getSumIte_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getSumIte_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getSumIte_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getSumIte_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getSumIte_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getSumIte_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getSumIte_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getSumIte_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ite_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Rec_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getSumRec_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getSumRec_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getSumRec_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getSumRec_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getSumRec_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getSumRec_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getSumRec_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getSumRec_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getSumRec_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getSumRec_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Rec_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FAB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getSumFAB_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getSumFAB_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getSumFAB_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getSumFAB_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getSumFAB_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getSumFAB_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getSumFAB_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getSumFAB_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getSumFAB_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getSumFAB_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FAB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NAB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getSumNAB_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getSumNAB_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getSumNAB_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getSumNAB_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getSumNAB_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getSumNAB_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getSumNAB_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getSumNAB_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getSumNAB_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getSumNAB_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NAB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define KAB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getSumKAB_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getSumKAB_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getSumKAB_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getSumKAB_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getSumKAB_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getSumKAB_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getSumKAB_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getSumKAB_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getSumKAB_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getSumKAB_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef KAB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getSum_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDot_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getDotDef_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getDotDef_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getDotDef_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getDotDef_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getDotDef_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDotDef_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDotDef_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDotDef_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDotDef_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDotDef_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ite_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getDotIte_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getDotIte_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getDotIte_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getDotIte_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getDotIte_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDotIte_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDotIte_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDotIte_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDotIte_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDotIte_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ite_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Rec_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getDotRec_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getDotRec_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getDotRec_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getDotRec_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getDotRec_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDotRec_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDotRec_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDotRec_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDotRec_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDotRec_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Rec_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FAB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getDotFAB_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getDotFAB_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getDotFAB_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getDotFAB_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getDotFAB_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDotFAB_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDotFAB_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDotFAB_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDotFAB_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDotFAB_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FAB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NAB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getDotNAB_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getDotNAB_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getDotNAB_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getDotNAB_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getDotNAB_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDotNAB_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDotNAB_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDotNAB_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDotNAB_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDotNAB_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NAB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define KAB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getDotKAB_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getDotKAB_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getDotKAB_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getDotKAB_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getDotKAB_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDotKAB_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDotKAB_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDotKAB_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDotKAB_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDotKAB_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef KAB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDot_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines