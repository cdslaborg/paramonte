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
!>  This file contains procedure implementations of [pm_mathCumSum](@ref pm_mathCumSum).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 PM, September 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_mathCumSum) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_arrayReverse, only: setReversed
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCumSum_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getCumSum_IK5
        use pm_kind, only: TKG => IK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getCumSum_IK4
        use pm_kind, only: TKG => IK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getCumSum_IK3
        use pm_kind, only: TKG => IK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getCumSum_IK2
        use pm_kind, only: TKG => IK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getCumSum_IK1
        use pm_kind, only: TKG => IK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getCumSum_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCumSum_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCumSum_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCumSum_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCumSum_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCumSum_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCumSum_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCumSum_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCumSum_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCumSum_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCumSum_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCumSum_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Old_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define For_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Non_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setCumSumOldDefDef_IK5
        use pm_kind, only: TKG => IK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setCumSumOldDefDef_IK4
        use pm_kind, only: TKG => IK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setCumSumOldDefDef_IK3
        use pm_kind, only: TKG => IK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setCumSumOldDefDef_IK2
        use pm_kind, only: TKG => IK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setCumSumOldDefDef_IK1
        use pm_kind, only: TKG => IK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCumSumOldDefDef_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCumSumOldDefDef_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCumSumOldDefDef_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCumSumOldDefDef_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCumSumOldDefDef_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCumSumOldDefDef_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCumSumOldDefDef_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCumSumOldDefDef_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCumSumOldDefDef_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCumSumOldDefDef_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Non_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef For_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Old_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define For_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Non_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setCumSumNewDefDef_IK5
        use pm_kind, only: TKG => IK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setCumSumNewDefDef_IK4
        use pm_kind, only: TKG => IK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setCumSumNewDefDef_IK3
        use pm_kind, only: TKG => IK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setCumSumNewDefDef_IK2
        use pm_kind, only: TKG => IK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setCumSumNewDefDef_IK1
        use pm_kind, only: TKG => IK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCumSumNewDefDef_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCumSumNewDefDef_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCumSumNewDefDef_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCumSumNewDefDef_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCumSumNewDefDef_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCumSumNewDefDef_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCumSumNewDefDef_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCumSumNewDefDef_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCumSumNewDefDef_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCumSumNewDefDef_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Non_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef For_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCumSum_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCumSum_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Old_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define For_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Non_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setCumSumOldForNon_IK5
        use pm_kind, only: TKG => IK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setCumSumOldForNon_IK4
        use pm_kind, only: TKG => IK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setCumSumOldForNon_IK3
        use pm_kind, only: TKG => IK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setCumSumOldForNon_IK2
        use pm_kind, only: TKG => IK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setCumSumOldForNon_IK1
        use pm_kind, only: TKG => IK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCumSumOldForNon_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCumSumOldForNon_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCumSumOldForNon_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCumSumOldForNon_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCumSumOldForNon_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCumSumOldForNon_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCumSumOldForNon_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCumSumOldForNon_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCumSumOldForNon_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCumSumOldForNon_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Non_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Rev_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setCumSumOldForRev_IK5
        use pm_kind, only: TKG => IK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setCumSumOldForRev_IK4
        use pm_kind, only: TKG => IK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setCumSumOldForRev_IK3
        use pm_kind, only: TKG => IK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setCumSumOldForRev_IK2
        use pm_kind, only: TKG => IK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setCumSumOldForRev_IK1
        use pm_kind, only: TKG => IK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCumSumOldForRev_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCumSumOldForRev_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCumSumOldForRev_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCumSumOldForRev_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCumSumOldForRev_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCumSumOldForRev_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCumSumOldForRev_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCumSumOldForRev_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCumSumOldForRev_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCumSumOldForRev_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Rev_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef For_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Bac_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Non_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setCumSumOldBacNon_IK5
        use pm_kind, only: TKG => IK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setCumSumOldBacNon_IK4
        use pm_kind, only: TKG => IK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setCumSumOldBacNon_IK3
        use pm_kind, only: TKG => IK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setCumSumOldBacNon_IK2
        use pm_kind, only: TKG => IK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setCumSumOldBacNon_IK1
        use pm_kind, only: TKG => IK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCumSumOldBacNon_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCumSumOldBacNon_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCumSumOldBacNon_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCumSumOldBacNon_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCumSumOldBacNon_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCumSumOldBacNon_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCumSumOldBacNon_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCumSumOldBacNon_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCumSumOldBacNon_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCumSumOldBacNon_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Non_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Rev_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setCumSumOldBacRev_IK5
        use pm_kind, only: TKG => IK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setCumSumOldBacRev_IK4
        use pm_kind, only: TKG => IK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setCumSumOldBacRev_IK3
        use pm_kind, only: TKG => IK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setCumSumOldBacRev_IK2
        use pm_kind, only: TKG => IK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setCumSumOldBacRev_IK1
        use pm_kind, only: TKG => IK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCumSumOldBacRev_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCumSumOldBacRev_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCumSumOldBacRev_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCumSumOldBacRev_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCumSumOldBacRev_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCumSumOldBacRev_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCumSumOldBacRev_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCumSumOldBacRev_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCumSumOldBacRev_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCumSumOldBacRev_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Rev_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Bac_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Old_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define For_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Non_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setCumSumNewForNon_IK5
        use pm_kind, only: TKG => IK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setCumSumNewForNon_IK4
        use pm_kind, only: TKG => IK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setCumSumNewForNon_IK3
        use pm_kind, only: TKG => IK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setCumSumNewForNon_IK2
        use pm_kind, only: TKG => IK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setCumSumNewForNon_IK1
        use pm_kind, only: TKG => IK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCumSumNewForNon_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCumSumNewForNon_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCumSumNewForNon_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCumSumNewForNon_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCumSumNewForNon_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCumSumNewForNon_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCumSumNewForNon_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCumSumNewForNon_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCumSumNewForNon_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCumSumNewForNon_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Non_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Rev_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setCumSumNewForRev_IK5
        use pm_kind, only: TKG => IK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setCumSumNewForRev_IK4
        use pm_kind, only: TKG => IK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setCumSumNewForRev_IK3
        use pm_kind, only: TKG => IK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setCumSumNewForRev_IK2
        use pm_kind, only: TKG => IK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setCumSumNewForRev_IK1
        use pm_kind, only: TKG => IK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCumSumNewForRev_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCumSumNewForRev_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCumSumNewForRev_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCumSumNewForRev_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCumSumNewForRev_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCumSumNewForRev_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCumSumNewForRev_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCumSumNewForRev_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCumSumNewForRev_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCumSumNewForRev_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Rev_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef For_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Bac_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Non_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setCumSumNewBacNon_IK5
        use pm_kind, only: TKG => IK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setCumSumNewBacNon_IK4
        use pm_kind, only: TKG => IK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setCumSumNewBacNon_IK3
        use pm_kind, only: TKG => IK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setCumSumNewBacNon_IK2
        use pm_kind, only: TKG => IK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setCumSumNewBacNon_IK1
        use pm_kind, only: TKG => IK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCumSumNewBacNon_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCumSumNewBacNon_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCumSumNewBacNon_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCumSumNewBacNon_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCumSumNewBacNon_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCumSumNewBacNon_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCumSumNewBacNon_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCumSumNewBacNon_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCumSumNewBacNon_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCumSumNewBacNon_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Non_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Rev_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setCumSumNewBacRev_IK5
        use pm_kind, only: TKG => IK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setCumSumNewBacRev_IK4
        use pm_kind, only: TKG => IK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setCumSumNewBacRev_IK3
        use pm_kind, only: TKG => IK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setCumSumNewBacRev_IK2
        use pm_kind, only: TKG => IK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setCumSumNewBacRev_IK1
        use pm_kind, only: TKG => IK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCumSumNewBacRev_CK5
        use pm_kind, only: TKG => CK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCumSumNewBacRev_CK4
        use pm_kind, only: TKG => CK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCumSumNewBacRev_CK3
        use pm_kind, only: TKG => CK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCumSumNewBacRev_CK2
        use pm_kind, only: TKG => CK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCumSumNewBacRev_CK1
        use pm_kind, only: TKG => CK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCumSumNewBacRev_RK5
        use pm_kind, only: TKG => RK5
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCumSumNewBacRev_RK4
        use pm_kind, only: TKG => RK4
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCumSumNewBacRev_RK3
        use pm_kind, only: TKG => RK3
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCumSumNewBacRev_RK2
        use pm_kind, only: TKG => RK2
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCumSumNewBacRev_RK1
        use pm_kind, only: TKG => RK1
#include "pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Rev_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Bac_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCumSum_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines
