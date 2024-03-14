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
!>  This file contains procedure implementations of [pm_complexMinMax](@ref pm_complexMinMax).
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_complexMinMax) routines ! LCOV_EXCL_LINE

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

#define min_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure min_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure min_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure min_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure min_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure min_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef min_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define max_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure max_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure max_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure max_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure max_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure max_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef max_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define minval_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure minvalALL_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure minvalALL_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure minvalALL_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure minvalALL_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure minvalALL_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure minvalALL_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure minvalALL_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure minvalALL_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure minvalALL_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure minvalALL_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DIM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure minvalDIM_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure minvalDIM_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure minvalDIM_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure minvalDIM_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure minvalDIM_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure minvalDIM_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure minvalDIM_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure minvalDIM_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure minvalDIM_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure minvalDIM_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DIM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef minval_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define maxval_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure maxvalALL_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure maxvalALL_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure maxvalALL_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure maxvalALL_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure maxvalALL_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure maxvalALL_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure maxvalALL_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure maxvalALL_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure maxvalALL_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure maxvalALL_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DIM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure maxvalDIM_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure maxvalDIM_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure maxvalDIM_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure maxvalDIM_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure maxvalDIM_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure maxvalDIM_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure maxvalDIM_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure maxvalDIM_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure maxvalDIM_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure maxvalDIM_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DIM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef maxval_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define minloc_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure minlocALL_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure minlocALL_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure minlocALL_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure minlocALL_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure minlocALL_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure minlocALL_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure minlocALL_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure minlocALL_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure minlocALL_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure minlocALL_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DIM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure minlocDIM_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure minlocDIM_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure minlocDIM_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure minlocDIM_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure minlocDIM_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure minlocDIM_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure minlocDIM_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure minlocDIM_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure minlocDIM_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure minlocDIM_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DIM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef minloc_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define maxloc_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure maxlocALL_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure maxlocALL_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure maxlocALL_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure maxlocALL_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure maxlocALL_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure maxlocALL_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure maxlocALL_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure maxlocALL_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure maxlocALL_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure maxlocALL_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DIM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure maxlocDIM_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure maxlocDIM_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure maxlocDIM_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure maxlocDIM_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure maxlocDIM_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure maxlocDIM_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure maxlocDIM_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure maxlocDIM_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure maxlocDIM_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure maxlocDIM_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_complexMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DIM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef maxloc_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  CHECK_ASSERTION

end submodule routines