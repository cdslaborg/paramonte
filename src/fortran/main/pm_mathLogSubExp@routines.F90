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
!>  This file contains procedure implementations of [pm_mathLogSubExp](@ref pm_mathLogSubExp).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_mathLogSubExp) routines ! LCOV_EXCL_LINE

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

#define getLogSubExp_ENABLED 1

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
    module procedure getLogSubExpSeqSL_CK5
        use pm_kind, only: CKC => CK5
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getLogSubExpSeqSL_CK4
        use pm_kind, only: CKC => CK4
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getLogSubExpSeqSL_CK3
        use pm_kind, only: CKC => CK3
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getLogSubExpSeqSL_CK2
        use pm_kind, only: CKC => CK2
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getLogSubExpSeqSL_CK1
        use pm_kind, only: CKC => CK1
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogSubExpSeqSL_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogSubExpSeqSL_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogSubExpSeqSL_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogSubExpSeqSL_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogSubExpSeqSL_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathLogSubExp@routines.inc.F90"
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
    module procedure getLogSubExpSelSL_CK5
        use pm_kind, only: CKC => CK5
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getLogSubExpSelSL_CK4
        use pm_kind, only: CKC => CK4
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getLogSubExpSelSL_CK3
        use pm_kind, only: CKC => CK3
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getLogSubExpSelSL_CK2
        use pm_kind, only: CKC => CK2
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getLogSubExpSelSL_CK1
        use pm_kind, only: CKC => CK1
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogSubExpSelSL_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogSubExpSelSL_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogSubExpSelSL_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogSubExpSelSL_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogSubExpSelSL_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathLogSubExp@routines.inc.F90"
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
    module procedure getLogSubExpSeqMM_CK5
        use pm_kind, only: CKC => CK5
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getLogSubExpSeqMM_CK4
        use pm_kind, only: CKC => CK4
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getLogSubExpSeqMM_CK3
        use pm_kind, only: CKC => CK3
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getLogSubExpSeqMM_CK2
        use pm_kind, only: CKC => CK2
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getLogSubExpSeqMM_CK1
        use pm_kind, only: CKC => CK1
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogSubExpSeqMM_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogSubExpSeqMM_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogSubExpSeqMM_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogSubExpSeqMM_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogSubExpSeqMM_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathLogSubExp@routines.inc.F90"
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
    module procedure getLogSubExpSelMM_CK5
        use pm_kind, only: CKC => CK5
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getLogSubExpSelMM_CK4
        use pm_kind, only: CKC => CK4
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getLogSubExpSelMM_CK3
        use pm_kind, only: CKC => CK3
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getLogSubExpSelMM_CK2
        use pm_kind, only: CKC => CK2
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getLogSubExpSelMM_CK1
        use pm_kind, only: CKC => CK1
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogSubExpSelMM_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogSubExpSelMM_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogSubExpSelMM_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogSubExpSelMM_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogSubExpSelMM_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathLogSubExp@routines.inc.F90"
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

#undef getLogSubExp_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines