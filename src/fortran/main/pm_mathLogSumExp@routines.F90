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
!>  This file contains procedure implementations of [pm_mathLogSumExp](@ref pm_mathLogSumExp).
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 PM, September 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_mathLogSumExp) routines ! LCOV_EXCL_LINE

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

#define getLogSumExp_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Seq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getLogSumExpDefSeq_CK5
        use pm_kind, only: TKC => CK5
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getLogSumExpDefSeq_CK4
        use pm_kind, only: TKC => CK4
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getLogSumExpDefSeq_CK3
        use pm_kind, only: TKC => CK3
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getLogSumExpDefSeq_CK2
        use pm_kind, only: TKC => CK2
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getLogSumExpDefSeq_CK1
        use pm_kind, only: TKC => CK1
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogSumExpDefSeq_RK5
        use pm_kind, only: TKC => RK5
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogSumExpDefSeq_RK4
        use pm_kind, only: TKC => RK4
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogSumExpDefSeq_RK3
        use pm_kind, only: TKC => RK3
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogSumExpDefSeq_RK2
        use pm_kind, only: TKC => RK2
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogSumExpDefSeq_RK1
        use pm_kind, only: TKC => RK1
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Seq_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Max_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Seq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getLogSumExpMaxSeq_CK5
        use pm_kind, only: TKC => CK5
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getLogSumExpMaxSeq_CK4
        use pm_kind, only: TKC => CK4
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getLogSumExpMaxSeq_CK3
        use pm_kind, only: TKC => CK3
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getLogSumExpMaxSeq_CK2
        use pm_kind, only: TKC => CK2
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getLogSumExpMaxSeq_CK1
        use pm_kind, only: TKC => CK1
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogSumExpMaxSeq_RK5
        use pm_kind, only: TKC => RK5
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogSumExpMaxSeq_RK4
        use pm_kind, only: TKC => RK4
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogSumExpMaxSeq_RK3
        use pm_kind, only: TKC => RK3
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogSumExpMaxSeq_RK2
        use pm_kind, only: TKC => RK2
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogSumExpMaxSeq_RK1
        use pm_kind, only: TKC => RK1
#include "pm_mathLogSumExp@routines.inc.F90"
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
    module procedure getLogSumExpMaxSel_CK5
        use pm_kind, only: TKC => CK5
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getLogSumExpMaxSel_CK4
        use pm_kind, only: TKC => CK4
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getLogSumExpMaxSel_CK3
        use pm_kind, only: TKC => CK3
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getLogSumExpMaxSel_CK2
        use pm_kind, only: TKC => CK2
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getLogSumExpMaxSel_CK1
        use pm_kind, only: TKC => CK1
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogSumExpMaxSel_RK5
        use pm_kind, only: TKC => RK5
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogSumExpMaxSel_RK4
        use pm_kind, only: TKC => RK4
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogSumExpMaxSel_RK3
        use pm_kind, only: TKC => RK3
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogSumExpMaxSel_RK2
        use pm_kind, only: TKC => RK2
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogSumExpMaxSel_RK1
        use pm_kind, only: TKC => RK1
#include "pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Sel_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Max_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogSumExp_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines