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
!>  This file contains procedure implementations of [pm_sampleCCF](@ref pm_sampleCCF).
!>
!>  final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_sampleCCF) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_fftpack, only: setFFTF, setFFTI, setFFTR
    use pm_arraySort, only: isAscending
    use pm_arrayResize, only: setResized
    use pm_sampleShift, only: setShifted
    use pm_sampleScale, only: setScaled
    use pm_sampleNorm, only: setNormed
    use pm_sampleMean, only: getMean

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getACF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getACF_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getACF_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getACF_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getACF_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getACF_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getACF_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getACF_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getACF_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getACF_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getACF_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getACF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setACF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setACF_FP_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setACF_FP_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setACF_FP_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setACF_FP_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setACF_FP_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setACF_FP_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setACF_FP_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setACF_FP_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setACF_FP_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setACF_FP_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setACF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCCF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FG_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getCCF_FG_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCCF_FG_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCCF_FG_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCCF_FG_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCCF_FG_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCCF_FG_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCCF_FG_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCCF_FG_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCCF_FG_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCCF_FG_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FG_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCCF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCCF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FG_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCCF_FP_FG_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCCF_FP_FG_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCCF_FP_FG_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCCF_FP_FG_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCCF_FP_FG_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCCF_FP_FG_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCCF_FP_FG_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCCF_FP_FG_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCCF_FP_FG_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCCF_FP_FG_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FG_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCCF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines