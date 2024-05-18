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
!>  This file contains the implementations of the tests of module [pm_arraySpace](@ref pm_arraySpace).
!>  
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Sunday 4:33 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_arraySpace) routines

    use pm_kind, only: LK
    use pm_val2str, only: getStr
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLinSpace_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getLinSpace_CK5
        use pm_kind, only: TKC => CK5
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getLinSpace_CK4
        use pm_kind, only: TKC => CK4
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getLinSpace_CK3
        use pm_kind, only: TKC => CK3
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getLinSpace_CK2
        use pm_kind, only: TKC => CK2
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getLinSpace_CK1
        use pm_kind, only: TKC => CK1
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getLinSpace_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getLinSpace_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getLinSpace_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getLinSpace_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getLinSpace_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLinSpace_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setLinSpace_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setLinSpace_CK5
        use pm_kind, only: TKC => CK5
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setLinSpace_CK4
        use pm_kind, only: TKC => CK4
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setLinSpace_CK3
        use pm_kind, only: TKC => CK3
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setLinSpace_CK2
        use pm_kind, only: TKC => CK2
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setLinSpace_CK1
        use pm_kind, only: TKC => CK1
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setLinSpace_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setLinSpace_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setLinSpace_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setLinSpace_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setLinSpace_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setLinSpace_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogSpace_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getLogSpace_CK5
        use pm_kind, only: TKC => CK5
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getLogSpace_CK4
        use pm_kind, only: TKC => CK4
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getLogSpace_CK3
        use pm_kind, only: TKC => CK3
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getLogSpace_CK2
        use pm_kind, only: TKC => CK2
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getLogSpace_CK1
        use pm_kind, only: TKC => CK1
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getLogSpace_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getLogSpace_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getLogSpace_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getLogSpace_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getLogSpace_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_arraySpace@routines.inc.F90"
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
    module procedure test_setLogSpace_CK5
        use pm_kind, only: TKC => CK5
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setLogSpace_CK4
        use pm_kind, only: TKC => CK4
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setLogSpace_CK3
        use pm_kind, only: TKC => CK3
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setLogSpace_CK2
        use pm_kind, only: TKC => CK2
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setLogSpace_CK1
        use pm_kind, only: TKC => CK1
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setLogSpace_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setLogSpace_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setLogSpace_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setLogSpace_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setLogSpace_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_arraySpace@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setLogSpace_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE