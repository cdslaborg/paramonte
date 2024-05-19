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
!>  This file contains procedure implementations of [pm_mathLogAddExp](@ref pm_mathLogAddExp).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_mathLogAddExp) routines ! LCOV_EXCL_LINE

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

#define getLogAddExp_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Seq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getLogAddExpSeqSL_CK5
        use pm_kind, only: CKG => CK5
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getLogAddExpSeqSL_CK4
        use pm_kind, only: CKG => CK4
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getLogAddExpSeqSL_CK3
        use pm_kind, only: CKG => CK3
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getLogAddExpSeqSL_CK2
        use pm_kind, only: CKG => CK2
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getLogAddExpSeqSL_CK1
        use pm_kind, only: CKG => CK1
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogAddExpSeqSL_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogAddExpSeqSL_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogAddExpSeqSL_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogAddExpSeqSL_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogAddExpSeqSL_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Seq_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Sel_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getLogAddExpSelSL_CK5
        use pm_kind, only: CKG => CK5
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getLogAddExpSelSL_CK4
        use pm_kind, only: CKG => CK4
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getLogAddExpSelSL_CK3
        use pm_kind, only: CKG => CK3
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getLogAddExpSelSL_CK2
        use pm_kind, only: CKG => CK2
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getLogAddExpSelSL_CK1
        use pm_kind, only: CKG => CK1
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogAddExpSelSL_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogAddExpSelSL_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogAddExpSelSL_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogAddExpSelSL_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogAddExpSelSL_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Sel_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Seq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getLogAddExpSeqMM_CK5
        use pm_kind, only: CKG => CK5
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getLogAddExpSeqMM_CK4
        use pm_kind, only: CKG => CK4
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getLogAddExpSeqMM_CK3
        use pm_kind, only: CKG => CK3
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getLogAddExpSeqMM_CK2
        use pm_kind, only: CKG => CK2
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getLogAddExpSeqMM_CK1
        use pm_kind, only: CKG => CK1
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogAddExpSeqMM_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogAddExpSeqMM_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogAddExpSeqMM_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogAddExpSeqMM_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogAddExpSeqMM_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Seq_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Sel_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getLogAddExpSelMM_CK5
        use pm_kind, only: CKG => CK5
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getLogAddExpSelMM_CK4
        use pm_kind, only: CKG => CK4
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getLogAddExpSelMM_CK3
        use pm_kind, only: CKG => CK3
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getLogAddExpSelMM_CK2
        use pm_kind, only: CKG => CK2
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getLogAddExpSelMM_CK1
        use pm_kind, only: CKG => CK1
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogAddExpSelMM_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogAddExpSelMM_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogAddExpSelMM_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogAddExpSelMM_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogAddExpSelMM_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Sel_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogAddExp_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines