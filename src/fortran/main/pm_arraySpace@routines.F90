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
!>  This file contains procedure implementations of [pm_arraySpace](@ref pm_arraySpace).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 PM, September 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_arraySpace) routines ! LCOV_EXCL_LINE

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

#define getLinSpace_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLinSpace_RK5
        use pm_kind, only: TKC => RK5
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLinSpace_RK4
        use pm_kind, only: TKC => RK4
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLinSpace_RK3
        use pm_kind, only: TKC => RK3
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLinSpace_RK2
        use pm_kind, only: TKC => RK2
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLinSpace_RK1
        use pm_kind, only: TKC => RK1
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getLinSpace_CK5
        use pm_kind, only: TKC => CK5
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getLinSpace_CK4
        use pm_kind, only: TKC => CK4
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getLinSpace_CK3
        use pm_kind, only: TKC => CK3
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getLinSpace_CK2
        use pm_kind, only: TKC => CK2
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getLinSpace_CK1
        use pm_kind, only: TKC => CK1
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLinSpace_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setLinSpace_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setLinSpace_RK5
        use pm_kind, only: TKC => RK5
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setLinSpace_RK4
        use pm_kind, only: TKC => RK4
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setLinSpace_RK3
        use pm_kind, only: TKC => RK3
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setLinSpace_RK2
        use pm_kind, only: TKC => RK2
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setLinSpace_RK1
        use pm_kind, only: TKC => RK1
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setLinSpace_CK5
        use pm_kind, only: TKC => CK5
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setLinSpace_CK4
        use pm_kind, only: TKC => CK4
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setLinSpace_CK3
        use pm_kind, only: TKC => CK3
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setLinSpace_CK2
        use pm_kind, only: TKC => CK2
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setLinSpace_CK1
        use pm_kind, only: TKC => CK1
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setLinSpace_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogSpace_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getLogSpace_CK5
        use pm_kind, only: TKC => CK5
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getLogSpace_CK4
        use pm_kind, only: TKC => CK4
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getLogSpace_CK3
        use pm_kind, only: TKC => CK3
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getLogSpace_CK2
        use pm_kind, only: TKC => CK2
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getLogSpace_CK1
        use pm_kind, only: TKC => CK1
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogSpace_RK5
        use pm_kind, only: TKC => RK5
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogSpace_RK4
        use pm_kind, only: TKC => RK4
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogSpace_RK3
        use pm_kind, only: TKC => RK3
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogSpace_RK2
        use pm_kind, only: TKC => RK2
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogSpace_RK1
        use pm_kind, only: TKC => RK1
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogSpace_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setLogSpace_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setLogSpace_CK5
        use pm_kind, only: TKC => CK5
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setLogSpace_CK4
        use pm_kind, only: TKC => CK4
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setLogSpace_CK3
        use pm_kind, only: TKC => CK3
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setLogSpace_CK2
        use pm_kind, only: TKC => CK2
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setLogSpace_CK1
        use pm_kind, only: TKC => CK1
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setLogSpace_RK5
        use pm_kind, only: TKC => RK5
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setLogSpace_RK4
        use pm_kind, only: TKC => RK4
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setLogSpace_RK3
        use pm_kind, only: TKC => RK3
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setLogSpace_RK2
        use pm_kind, only: TKC => RK2
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setLogSpace_RK1
        use pm_kind, only: TKC => RK1
#include "pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setLogSpace_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines