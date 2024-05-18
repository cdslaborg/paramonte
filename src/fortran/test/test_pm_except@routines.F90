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

!>  \brief This file contains the implementations of the tests of module [pm_except](@ref pm_except).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_except) routines

    use pm_distUnif, only: setUnifRand
    use pm_val2str, only: getStr
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getInfPos_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getInfPos_CK5_1
        use pm_kind, only: CKC => CK5
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getInfPos_CK4_1
        use pm_kind, only: CKC => CK4
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getInfPos_CK3_1
        use pm_kind, only: CKC => CK3
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getInfPos_CK2_1
        use pm_kind, only: CKC => CK2
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getInfPos_CK1_1
        use pm_kind, only: CKC => CK1
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getInfPos_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getInfPos_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getInfPos_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getInfPos_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getInfPos_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getInfPos_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setInfPos_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setInfPos_CK5_1
        use pm_kind, only: CKC => CK5
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setInfPos_CK4_1
        use pm_kind, only: CKC => CK4
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setInfPos_CK3_1
        use pm_kind, only: CKC => CK3
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setInfPos_CK2_1
        use pm_kind, only: CKC => CK2
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setInfPos_CK1_1
        use pm_kind, only: CKC => CK1
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setInfPos_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setInfPos_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setInfPos_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setInfPos_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setInfPos_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setInfPos_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getInfNeg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getInfNeg_CK5_1
        use pm_kind, only: CKC => CK5
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getInfNeg_CK4_1
        use pm_kind, only: CKC => CK4
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getInfNeg_CK3_1
        use pm_kind, only: CKC => CK3
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getInfNeg_CK2_1
        use pm_kind, only: CKC => CK2
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getInfNeg_CK1_1
        use pm_kind, only: CKC => CK1
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getInfNeg_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getInfNeg_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getInfNeg_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getInfNeg_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getInfNeg_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getInfNeg_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setInfNeg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setInfNeg_CK5_1
        use pm_kind, only: CKC => CK5
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setInfNeg_CK4_1
        use pm_kind, only: CKC => CK4
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setInfNeg_CK3_1
        use pm_kind, only: CKC => CK3
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setInfNeg_CK2_1
        use pm_kind, only: CKC => CK2
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setInfNeg_CK1_1
        use pm_kind, only: CKC => CK1
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setInfNeg_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setInfNeg_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setInfNeg_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setInfNeg_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setInfNeg_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setInfNeg_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNAN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getNAN_CK5_1
        use pm_kind, only: CKC => CK5
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getNAN_CK4_1
        use pm_kind, only: CKC => CK4
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getNAN_CK3_1
        use pm_kind, only: CKC => CK3
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getNAN_CK2_1
        use pm_kind, only: CKC => CK2
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getNAN_CK1_1
        use pm_kind, only: CKC => CK1
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getNAN_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getNAN_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getNAN_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getNAN_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getNAN_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setNAN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setNAN_CK5_1
        use pm_kind, only: CKC => CK5
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setNAN_CK4_1
        use pm_kind, only: CKC => CK4
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setNAN_CK3_1
        use pm_kind, only: CKC => CK3
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setNAN_CK2_1
        use pm_kind, only: CKC => CK2
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setNAN_CK1_1
        use pm_kind, only: CKC => CK1
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setNAN_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setNAN_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setNAN_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setNAN_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setNAN_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_except@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setNAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE