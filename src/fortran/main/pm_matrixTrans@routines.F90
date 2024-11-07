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
!>  This file contains procedure implementations of [pm_matrixTrans](@ref pm_matrixTrans).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_matrixTrans) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_matrixLUP, only: setMatLUP
   !use pm_matrixInit, only: setMatInit
    use pm_matrixCopy, only: setMatCopy, rdpack, lfpack, dia
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatTrans_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Fix_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Symm_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Old_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatTransSymmOldFix_SK5
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatTransSymmOldFix_SK4
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatTransSymmOldFix_SK3
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatTransSymmOldFix_SK2
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatTransSymmOldFix_SK1
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatTransSymmOldFix_IK5
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatTransSymmOldFix_IK4
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatTransSymmOldFix_IK3
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatTransSymmOldFix_IK2
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatTransSymmOldFix_IK1
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatTransSymmOldFix_LK5
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatTransSymmOldFix_LK4
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatTransSymmOldFix_LK3
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatTransSymmOldFix_LK2
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatTransSymmOldFix_LK1
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatTransSymmOldFix_CK5
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatTransSymmOldFix_CK4
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatTransSymmOldFix_CK3
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatTransSymmOldFix_CK2
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatTransSymmOldFix_CK1
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatTransSymmOldFix_RK5
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatTransSymmOldFix_RK4
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatTransSymmOldFix_RK3
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatTransSymmOldFix_RK2
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatTransSymmOldFix_RK1
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Old_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatTransSymmNewFix_SK5
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatTransSymmNewFix_SK4
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatTransSymmNewFix_SK3
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatTransSymmNewFix_SK2
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatTransSymmNewFix_SK1
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatTransSymmNewFix_IK5
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatTransSymmNewFix_IK4
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatTransSymmNewFix_IK3
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatTransSymmNewFix_IK2
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatTransSymmNewFix_IK1
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatTransSymmNewFix_LK5
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatTransSymmNewFix_LK4
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatTransSymmNewFix_LK3
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatTransSymmNewFix_LK2
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatTransSymmNewFix_LK1
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatTransSymmNewFix_CK5
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatTransSymmNewFix_CK4
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatTransSymmNewFix_CK3
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatTransSymmNewFix_CK2
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatTransSymmNewFix_CK1
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatTransSymmNewFix_RK5
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatTransSymmNewFix_RK4
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatTransSymmNewFix_RK3
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatTransSymmNewFix_RK2
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatTransSymmNewFix_RK1
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Symm_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Herm_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Old_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatTransHermOldFix_CK5
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatTransHermOldFix_CK4
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatTransHermOldFix_CK3
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatTransHermOldFix_CK2
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatTransHermOldFix_CK1
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Old_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatTransHermNewFix_CK5
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatTransHermNewFix_CK4
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatTransHermNewFix_CK3
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatTransHermNewFix_CK2
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatTransHermNewFix_CK1
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Herm_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Fix_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arb_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Symm_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Old_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatTransSymmOldArb_SK5
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatTransSymmOldArb_SK4
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatTransSymmOldArb_SK3
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatTransSymmOldArb_SK2
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatTransSymmOldArb_SK1
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatTransSymmOldArb_IK5
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatTransSymmOldArb_IK4
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatTransSymmOldArb_IK3
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatTransSymmOldArb_IK2
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatTransSymmOldArb_IK1
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatTransSymmOldArb_LK5
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatTransSymmOldArb_LK4
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatTransSymmOldArb_LK3
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatTransSymmOldArb_LK2
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatTransSymmOldArb_LK1
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatTransSymmOldArb_CK5
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatTransSymmOldArb_CK4
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatTransSymmOldArb_CK3
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatTransSymmOldArb_CK2
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatTransSymmOldArb_CK1
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatTransSymmOldArb_RK5
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatTransSymmOldArb_RK4
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatTransSymmOldArb_RK3
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatTransSymmOldArb_RK2
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatTransSymmOldArb_RK1
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Old_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMatTransSymmNewArb_SK5
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMatTransSymmNewArb_SK4
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMatTransSymmNewArb_SK3
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMatTransSymmNewArb_SK2
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMatTransSymmNewArb_SK1
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatTransSymmNewArb_IK5
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatTransSymmNewArb_IK4
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatTransSymmNewArb_IK3
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatTransSymmNewArb_IK2
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatTransSymmNewArb_IK1
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMatTransSymmNewArb_LK5
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMatTransSymmNewArb_LK4
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMatTransSymmNewArb_LK3
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMatTransSymmNewArb_LK2
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMatTransSymmNewArb_LK1
        use pm_logicalCompare, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatTransSymmNewArb_CK5
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatTransSymmNewArb_CK4
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatTransSymmNewArb_CK3
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatTransSymmNewArb_CK2
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatTransSymmNewArb_CK1
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatTransSymmNewArb_RK5
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatTransSymmNewArb_RK4
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatTransSymmNewArb_RK3
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatTransSymmNewArb_RK2
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatTransSymmNewArb_RK1
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Symm_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Herm_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Old_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatTransHermOldArb_CK5
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatTransHermOldArb_CK4
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatTransHermOldArb_CK3
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatTransHermOldArb_CK2
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatTransHermOldArb_CK1
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Old_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatTransHermNewArb_CK5
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatTransHermNewArb_CK4
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatTransHermNewArb_CK3
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatTransHermNewArb_CK2
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatTransHermNewArb_CK1
        use pm_complexCompareLex, only: operator(<)
#include "pm_matrixTrans@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Herm_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arb_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatTrans_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines