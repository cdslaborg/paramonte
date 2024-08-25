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
!>  This file contains procedure implementations of [pm_matrixChol](@ref pm_matrixChol).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_matrixChol) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif
    use pm_matrixCopy, only: setMatCopy
    use pm_matrixInit, only: setMatInit
    use pm_matrixMulAdd, only: setMatMulAdd
    use pm_matrixSubset, only: uppDia, lowDia
    use pm_matrixClass, only: upperDiag, lowerDiag
    use pm_matrixTrans, only: transSymm, transHerm, transUnit
    use pm_matrixUpdate, only: setMatUpdate, hermitian, trans
    use pm_matrixMulTri, only: setMatMulTri
    !use pm_mathSum, only: dot_product => getDot

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setChoLow_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setChoLow_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setChoLow_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setChoLow_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setChoLow_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setChoLow_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setChoLow_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatChol_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatChol_UXD_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatChol_UXD_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatChol_UXD_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatChol_UXD_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatChol_UXD_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatChol_UXD_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatChol_UXD_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatChol_UXD_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatChol_UXD_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatChol_UXD_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatChol_XLD_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatChol_XLD_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatChol_XLD_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatChol_XLD_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatChol_XLD_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatChol_XLD_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatChol_XLD_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatChol_XLD_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatChol_XLD_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatChol_XLD_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMatChol_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatChol_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IMP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_IMP_AXX_UXD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_IMP_AXX_UXD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_IMP_AXX_UXD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_IMP_AXX_UXD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_IMP_AXX_UXD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_IMP_AXX_UXD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_IMP_AXX_UXD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_IMP_AXX_UXD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_IMP_AXX_UXD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_IMP_AXX_UXD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_IMP_AXX_XLD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_IMP_AXX_XLD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_IMP_AXX_XLD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_IMP_AXX_XLD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_IMP_AXX_XLD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_IMP_AXX_XLD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_IMP_AXX_XLD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_IMP_AXX_XLD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_IMP_AXX_XLD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_IMP_AXX_XLD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ABI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_IMP_ABI_UXD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_IMP_ABI_UXD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_IMP_ABI_UXD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_IMP_ABI_UXD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_IMP_ABI_UXD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_IMP_ABI_UXD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_IMP_ABI_UXD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_IMP_ABI_UXD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_IMP_ABI_UXD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_IMP_ABI_UXD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_IMP_ABI_XLD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_IMP_ABI_XLD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_IMP_ABI_XLD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_IMP_ABI_XLD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_IMP_ABI_XLD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_IMP_ABI_XLD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_IMP_ABI_XLD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_IMP_ABI_XLD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_IMP_ABI_XLD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_IMP_ABI_XLD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ABI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ABR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_IMP_ABR_UXD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_IMP_ABR_UXD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_IMP_ABR_UXD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_IMP_ABR_UXD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_IMP_ABR_UXD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_IMP_ABR_UXD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_IMP_ABR_UXD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_IMP_ABR_UXD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_IMP_ABR_UXD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_IMP_ABR_UXD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_IMP_ABR_XLD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_IMP_ABR_XLD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_IMP_ABR_XLD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_IMP_ABR_XLD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_IMP_ABR_XLD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_IMP_ABR_XLD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_IMP_ABR_XLD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_IMP_ABR_XLD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_IMP_ABR_XLD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_IMP_ABR_XLD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ABR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IMP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatChol_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatChol_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IMP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ANI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_IMP_ANI_UXD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_IMP_ANI_UXD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_IMP_ANI_UXD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_IMP_ANI_UXD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_IMP_ANI_UXD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_IMP_ANI_UXD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_IMP_ANI_UXD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_IMP_ANI_UXD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_IMP_ANI_UXD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_IMP_ANI_UXD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_IMP_ANI_XLD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_IMP_ANI_XLD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_IMP_ANI_XLD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_IMP_ANI_XLD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_IMP_ANI_XLD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_IMP_ANI_XLD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_IMP_ANI_XLD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_IMP_ANI_XLD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_IMP_ANI_XLD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_IMP_ANI_XLD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTH_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_IMP_ANI_UXD_OTH_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_IMP_ANI_UXD_OTH_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_IMP_ANI_UXD_OTH_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_IMP_ANI_UXD_OTH_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_IMP_ANI_UXD_OTH_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_IMP_ANI_UXD_OTH_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_IMP_ANI_UXD_OTH_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_IMP_ANI_UXD_OTH_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_IMP_ANI_UXD_OTH_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_IMP_ANI_UXD_OTH_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_IMP_ANI_XLD_OTH_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_IMP_ANI_XLD_OTH_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_IMP_ANI_XLD_OTH_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_IMP_ANI_XLD_OTH_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_IMP_ANI_XLD_OTH_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_IMP_ANI_XLD_OTH_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_IMP_ANI_XLD_OTH_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_IMP_ANI_XLD_OTH_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_IMP_ANI_XLD_OTH_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_IMP_ANI_XLD_OTH_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTH_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ANI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IMP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatChol_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatChol_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define EXP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_EXP_AXX_UXD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_EXP_AXX_UXD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_EXP_AXX_UXD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_EXP_AXX_UXD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_EXP_AXX_UXD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_EXP_AXX_UXD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_EXP_AXX_UXD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_EXP_AXX_UXD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_EXP_AXX_UXD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_EXP_AXX_UXD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_EXP_AXX_XLD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_EXP_AXX_XLD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_EXP_AXX_XLD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_EXP_AXX_XLD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_EXP_AXX_XLD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_EXP_AXX_XLD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_EXP_AXX_XLD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_EXP_AXX_XLD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_EXP_AXX_XLD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_EXP_AXX_XLD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ABI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_EXP_ABI_UXD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_EXP_ABI_UXD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_EXP_ABI_UXD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_EXP_ABI_UXD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_EXP_ABI_UXD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_EXP_ABI_UXD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_EXP_ABI_UXD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_EXP_ABI_UXD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_EXP_ABI_UXD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_EXP_ABI_UXD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_EXP_ABI_XLD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_EXP_ABI_XLD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_EXP_ABI_XLD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_EXP_ABI_XLD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_EXP_ABI_XLD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_EXP_ABI_XLD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_EXP_ABI_XLD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_EXP_ABI_XLD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_EXP_ABI_XLD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_EXP_ABI_XLD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ABI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ABR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_EXP_ABR_UXD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_EXP_ABR_UXD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_EXP_ABR_UXD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_EXP_ABR_UXD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_EXP_ABR_UXD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_EXP_ABR_UXD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_EXP_ABR_UXD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_EXP_ABR_UXD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_EXP_ABR_UXD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_EXP_ABR_UXD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_EXP_ABR_XLD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_EXP_ABR_XLD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_EXP_ABR_XLD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_EXP_ABR_XLD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_EXP_ABR_XLD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_EXP_ABR_XLD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_EXP_ABR_XLD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_EXP_ABR_XLD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_EXP_ABR_XLD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_EXP_ABR_XLD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ABR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef EXP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatChol_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatChol_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define EXP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ANI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_EXP_ANI_UXD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_EXP_ANI_UXD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_EXP_ANI_UXD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_EXP_ANI_UXD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_EXP_ANI_UXD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_EXP_ANI_UXD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_EXP_ANI_UXD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_EXP_ANI_UXD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_EXP_ANI_UXD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_EXP_ANI_UXD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_EXP_ANI_XLD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_EXP_ANI_XLD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_EXP_ANI_XLD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_EXP_ANI_XLD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_EXP_ANI_XLD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_EXP_ANI_XLD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_EXP_ANI_XLD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_EXP_ANI_XLD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_EXP_ANI_XLD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_EXP_ANI_XLD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTH_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_EXP_ANI_UXD_OTH_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_EXP_ANI_UXD_OTH_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_EXP_ANI_UXD_OTH_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_EXP_ANI_UXD_OTH_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_EXP_ANI_UXD_OTH_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_EXP_ANI_UXD_OTH_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_EXP_ANI_UXD_OTH_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_EXP_ANI_UXD_OTH_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_EXP_ANI_UXD_OTH_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_EXP_ANI_UXD_OTH_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMC_EXP_ANI_XLD_OTH_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMC_EXP_ANI_XLD_OTH_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMC_EXP_ANI_XLD_OTH_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMC_EXP_ANI_XLD_OTH_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMC_EXP_ANI_XLD_OTH_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMC_EXP_ANI_XLD_OTH_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMC_EXP_ANI_XLD_OTH_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMC_EXP_ANI_XLD_OTH_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMC_EXP_ANI_XLD_OTH_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMC_EXP_ANI_XLD_OTH_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixChol@routines.inc.F90"
#undef DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTH_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ANI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef EXP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatChol_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
