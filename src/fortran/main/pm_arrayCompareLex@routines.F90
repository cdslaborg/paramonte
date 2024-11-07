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
!>  This file contains procedure implementations of [pm_arrayCompareLex](@ref pm_arrayCompareLex).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_arrayCompareLex) routines ! LCOV_EXCL_LINE

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

#define isllt_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isllt_D0_D0_SK5
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isllt_D0_D0_SK4
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isllt_D0_D0_SK3
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isllt_D0_D0_SK2
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isllt_D0_D0_SK1
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isllt_D1_D1_SK5
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isllt_D1_D1_SK4
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isllt_D1_D1_SK3
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isllt_D1_D1_SK2
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isllt_D1_D1_SK1
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure isllt_D1_D1_IK5
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure isllt_D1_D1_IK4
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure isllt_D1_D1_IK3
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure isllt_D1_D1_IK2
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure isllt_D1_D1_IK1
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure isllt_D1_D1_LK5
        use pm_logicalCompare, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure isllt_D1_D1_LK4
        use pm_logicalCompare, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure isllt_D1_D1_LK3
        use pm_logicalCompare, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure isllt_D1_D1_LK2
        use pm_logicalCompare, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure isllt_D1_D1_LK1
        use pm_logicalCompare, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isllt_D1_D1_CK5
        use pm_complexCompareLex, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isllt_D1_D1_CK4
        use pm_complexCompareLex, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isllt_D1_D1_CK3
        use pm_complexCompareLex, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isllt_D1_D1_CK2
        use pm_complexCompareLex, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isllt_D1_D1_CK1
        use pm_complexCompareLex, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isllt_D1_D1_RK5
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isllt_D1_D1_RK4
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isllt_D1_D1_RK3
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isllt_D1_D1_RK2
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isllt_D1_D1_RK1
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isllt_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define islle_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure islle_D0_D0_SK5
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure islle_D0_D0_SK4
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure islle_D0_D0_SK3
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure islle_D0_D0_SK2
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure islle_D0_D0_SK1
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure islle_D1_D1_SK5
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure islle_D1_D1_SK4
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure islle_D1_D1_SK3
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure islle_D1_D1_SK2
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure islle_D1_D1_SK1
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure islle_D1_D1_IK5
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure islle_D1_D1_IK4
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure islle_D1_D1_IK3
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure islle_D1_D1_IK2
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure islle_D1_D1_IK1
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure islle_D1_D1_LK5
        use pm_logicalCompare, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure islle_D1_D1_LK4
        use pm_logicalCompare, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure islle_D1_D1_LK3
        use pm_logicalCompare, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure islle_D1_D1_LK2
        use pm_logicalCompare, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure islle_D1_D1_LK1
        use pm_logicalCompare, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure islle_D1_D1_CK5
        use pm_complexCompareLex, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure islle_D1_D1_CK4
        use pm_complexCompareLex, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure islle_D1_D1_CK3
        use pm_complexCompareLex, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure islle_D1_D1_CK2
        use pm_complexCompareLex, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure islle_D1_D1_CK1
        use pm_complexCompareLex, only: operator(<)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure islle_D1_D1_RK5
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure islle_D1_D1_RK4
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure islle_D1_D1_RK3
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure islle_D1_D1_RK2
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure islle_D1_D1_RK1
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef islle_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define islge_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure islge_D0_D0_SK5
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure islge_D0_D0_SK4
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure islge_D0_D0_SK3
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure islge_D0_D0_SK2
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure islge_D0_D0_SK1
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure islge_D1_D1_SK5
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure islge_D1_D1_SK4
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure islge_D1_D1_SK3
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure islge_D1_D1_SK2
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure islge_D1_D1_SK1
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure islge_D1_D1_IK5
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure islge_D1_D1_IK4
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure islge_D1_D1_IK3
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure islge_D1_D1_IK2
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure islge_D1_D1_IK1
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure islge_D1_D1_LK5
        use pm_logicalCompare, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure islge_D1_D1_LK4
        use pm_logicalCompare, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure islge_D1_D1_LK3
        use pm_logicalCompare, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure islge_D1_D1_LK2
        use pm_logicalCompare, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure islge_D1_D1_LK1
        use pm_logicalCompare, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure islge_D1_D1_CK5
        use pm_complexCompareLex, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure islge_D1_D1_CK4
        use pm_complexCompareLex, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure islge_D1_D1_CK3
        use pm_complexCompareLex, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure islge_D1_D1_CK2
        use pm_complexCompareLex, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure islge_D1_D1_CK1
        use pm_complexCompareLex, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure islge_D1_D1_RK5
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure islge_D1_D1_RK4
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure islge_D1_D1_RK3
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure islge_D1_D1_RK2
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure islge_D1_D1_RK1
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef islge_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define islgt_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure islgt_D0_D0_SK5
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure islgt_D0_D0_SK4
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure islgt_D0_D0_SK3
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure islgt_D0_D0_SK2
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure islgt_D0_D0_SK1
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure islgt_D1_D1_SK5
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure islgt_D1_D1_SK4
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure islgt_D1_D1_SK3
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure islgt_D1_D1_SK2
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure islgt_D1_D1_SK1
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure islgt_D1_D1_IK5
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure islgt_D1_D1_IK4
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure islgt_D1_D1_IK3
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure islgt_D1_D1_IK2
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure islgt_D1_D1_IK1
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure islgt_D1_D1_LK5
        use pm_logicalCompare, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure islgt_D1_D1_LK4
        use pm_logicalCompare, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure islgt_D1_D1_LK3
        use pm_logicalCompare, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure islgt_D1_D1_LK2
        use pm_logicalCompare, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure islgt_D1_D1_LK1
        use pm_logicalCompare, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure islgt_D1_D1_CK5
        use pm_complexCompareLex, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure islgt_D1_D1_CK4
        use pm_complexCompareLex, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure islgt_D1_D1_CK3
        use pm_complexCompareLex, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure islgt_D1_D1_CK2
        use pm_complexCompareLex, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure islgt_D1_D1_CK1
        use pm_complexCompareLex, only: operator(>)
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure islgt_D1_D1_RK5
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure islgt_D1_D1_RK4
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure islgt_D1_D1_RK3
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure islgt_D1_D1_RK2
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure islgt_D1_D1_RK1
#include "pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef islgt_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    module procedure islct_UP
!        use pm_kind, only: IK, LK
!        integer(IK) :: i
!        do i = 1_IK, min(GET_SIZE(array1, kind = IK), GET_SIZE(array2, kind = IK))
!            if (array1(i) EQUIVALENT_TO array2(i)) cycle
!            compares = array1(i) COMPARABLE_TO array2(i)
!            return
!        end do
!        compares = logical(GET_SIZE(array1, kind = IK) & ! LCOV_EXCL_LINE
!#if     isllt_ENABLED
!         < & ! LCOV_EXCL_LINE
!#elif   islle_ENABLED
!         <= & ! LCOV_EXCL_LINE
!#elif   islge_ENABLED
!         >= & ! LCOV_EXCL_LINE
!#elif   islgt_ENABLED
!         > & ! LCOV_EXCL_LINE
!#endif
!        GET_SIZE(array2, kind = IK), kind = LK)
!    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines