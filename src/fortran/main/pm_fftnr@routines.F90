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
!>  This file contains procedure implementations of [pm_fftnr](@ref pm_fftnr).
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 PM, September 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_fftnr) routines

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_mathExp, only: isIntPow
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getFFTF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getFFTF_CK5
        use pm_kind, only: TKC => CK5
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getFFTF_CK4
        use pm_kind, only: TKC => CK4
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getFFTF_CK3
        use pm_kind, only: TKC => CK3
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getFFTF_CK2
        use pm_kind, only: TKC => CK2
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getFFTF_CK1
        use pm_kind, only: TKC => CK1
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getFFTF_RK5
        use pm_kind, only: TKC => RK5
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getFFTF_RK4
        use pm_kind, only: TKC => RK4
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getFFTF_RK3
        use pm_kind, only: TKC => RK3
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getFFTF_RK2
        use pm_kind, only: TKC => RK2
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getFFTF_RK1
        use pm_kind, only: TKC => RK1
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getFFTF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getFFTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getFFTR_CK5
        use pm_kind, only: TKC => CK5
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getFFTR_CK4
        use pm_kind, only: TKC => CK4
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getFFTR_CK3
        use pm_kind, only: TKC => CK3
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getFFTR_CK2
        use pm_kind, only: TKC => CK2
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getFFTR_CK1
        use pm_kind, only: TKC => CK1
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getFFTR_RK5
        use pm_kind, only: TKC => RK5
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getFFTR_RK4
        use pm_kind, only: TKC => RK4
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getFFTR_RK3
        use pm_kind, only: TKC => RK3
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getFFTR_RK2
        use pm_kind, only: TKC => RK2
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getFFTR_RK1
        use pm_kind, only: TKC => RK1
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getFFTR_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getFFTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getFFTI_CK5
        use pm_kind, only: TKC => CK5
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getFFTI_CK4
        use pm_kind, only: TKC => CK4
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getFFTI_CK3
        use pm_kind, only: TKC => CK3
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getFFTI_CK2
        use pm_kind, only: TKC => CK2
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getFFTI_CK1
        use pm_kind, only: TKC => CK1
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getFFTI_RK5
        use pm_kind, only: TKC => RK5
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getFFTI_RK4
        use pm_kind, only: TKC => RK4
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getFFTI_RK3
        use pm_kind, only: TKC => RK3
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getFFTI_RK2
        use pm_kind, only: TKC => RK2
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getFFTI_RK1
        use pm_kind, only: TKC => RK1
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getFFTI_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setFFTF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setFFTF_CK5
        use pm_kind, only: TKC => CK5
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setFFTF_CK4
        use pm_kind, only: TKC => CK4
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setFFTF_CK3
        use pm_kind, only: TKC => CK3
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setFFTF_CK2
        use pm_kind, only: TKC => CK2
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setFFTF_CK1
        use pm_kind, only: TKC => CK1
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setFFTF_RK5
        use pm_kind, only: TKC => RK5
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setFFTF_RK4
        use pm_kind, only: TKC => RK4
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setFFTF_RK3
        use pm_kind, only: TKC => RK3
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setFFTF_RK2
        use pm_kind, only: TKC => RK2
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setFFTF_RK1
        use pm_kind, only: TKC => RK1
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setFFTF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setFFTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setFFTR_CK5
        use pm_kind, only: TKC => CK5
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setFFTR_CK4
        use pm_kind, only: TKC => CK4
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setFFTR_CK3
        use pm_kind, only: TKC => CK3
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setFFTR_CK2
        use pm_kind, only: TKC => CK2
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setFFTR_CK1
        use pm_kind, only: TKC => CK1
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setFFTR_RK5
        use pm_kind, only: TKC => RK5
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setFFTR_RK4
        use pm_kind, only: TKC => RK4
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setFFTR_RK3
        use pm_kind, only: TKC => RK3
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setFFTR_RK2
        use pm_kind, only: TKC => RK2
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setFFTR_RK1
        use pm_kind, only: TKC => RK1
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setFFTR_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setFFTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setFFTI_CK5
        use pm_kind, only: TKC => CK5
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setFFTI_CK4
        use pm_kind, only: TKC => CK4
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setFFTI_CK3
        use pm_kind, only: TKC => CK3
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setFFTI_CK2
        use pm_kind, only: TKC => CK2
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setFFTI_CK1
        use pm_kind, only: TKC => CK1
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setFFTI_RK5
        use pm_kind, only: TKC => RK5
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setFFTI_RK4
        use pm_kind, only: TKC => RK4
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setFFTI_RK3
        use pm_kind, only: TKC => RK3
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setFFTI_RK2
        use pm_kind, only: TKC => RK2
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setFFTI_RK1
        use pm_kind, only: TKC => RK1
#include "pm_fftnr@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setFFTI_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines