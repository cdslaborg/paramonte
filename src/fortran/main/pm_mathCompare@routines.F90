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
!>  This file contains procedure implementations of [pm_mathCompare](@ref pm_mathCompare).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 PM, September 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_mathCompare) routines ! LCOV_EXCL_LINE

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

#define isClose_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCloseDefault_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCloseDefault_CK_ENABLED 1

#if CK5_ENABLED
    module procedure isCloseDefault_CK5
        use pm_kind, only: CKG => CK5
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isCloseDefault_CK4
        use pm_kind, only: CKG => CK4
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isCloseDefault_CK3
        use pm_kind, only: CKG => CK3
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isCloseDefault_CK2
        use pm_kind, only: CKG => CK2
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isCloseDefault_CK1
        use pm_kind, only: CKG => CK1
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#undef isCloseDefault_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCloseDefault_RK_ENABLED 1

#if RK5_ENABLED
    module procedure isCloseDefault_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isCloseDefault_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isCloseDefault_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isCloseDefault_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isCloseDefault_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#undef isCloseDefault_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isCloseDefault_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCloseReference_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCloseReference_CK_ENABLED 1

#if CK5_ENABLED
    module procedure isCloseReference_CK5
        use pm_kind, only: CKG => CK5
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isCloseReference_CK4
        use pm_kind, only: CKG => CK4
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isCloseReference_CK3
        use pm_kind, only: CKG => CK3
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isCloseReference_CK2
        use pm_kind, only: CKG => CK2
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isCloseReference_CK1
        use pm_kind, only: CKG => CK1
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#undef isCloseReference_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCloseReference_RK_ENABLED 1

#if RK5_ENABLED
    module procedure isCloseReference_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isCloseReference_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isCloseReference_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isCloseReference_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isCloseReference_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#undef isCloseReference_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isCloseReference_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCloseStrong_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCloseStrong_CK_ENABLED 1

#if CK5_ENABLED
    module procedure isCloseStrong_CK5
        use pm_kind, only: CKG => CK5
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isCloseStrong_CK4
        use pm_kind, only: CKG => CK4
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isCloseStrong_CK3
        use pm_kind, only: CKG => CK3
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isCloseStrong_CK2
        use pm_kind, only: CKG => CK2
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isCloseStrong_CK1
        use pm_kind, only: CKG => CK1
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#undef isCloseStrong_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCloseStrong_RK_ENABLED 1

#if RK5_ENABLED
    module procedure isCloseStrong_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isCloseStrong_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isCloseStrong_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isCloseStrong_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isCloseStrong_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#undef isCloseStrong_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isCloseStrong_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCloseWeak_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCloseWeak_CK_ENABLED 1

#if CK5_ENABLED
    module procedure isCloseWeak_CK5
        use pm_kind, only: CKG => CK5
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isCloseWeak_CK4
        use pm_kind, only: CKG => CK4
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isCloseWeak_CK3
        use pm_kind, only: CKG => CK3
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isCloseWeak_CK2
        use pm_kind, only: CKG => CK2
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isCloseWeak_CK1
        use pm_kind, only: CKG => CK1
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#undef isCloseWeak_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCloseWeak_RK_ENABLED 1

#if RK5_ENABLED
    module procedure isCloseWeak_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isCloseWeak_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isCloseWeak_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isCloseWeak_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isCloseWeak_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#undef isCloseWeak_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isCloseWeak_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCloseMean_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCloseMean_CK_ENABLED 1

#if CK5_ENABLED
    module procedure isCloseMean_CK5
        use pm_kind, only: CKG => CK5
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isCloseMean_CK4
        use pm_kind, only: CKG => CK4
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isCloseMean_CK3
        use pm_kind, only: CKG => CK3
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isCloseMean_CK2
        use pm_kind, only: CKG => CK2
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isCloseMean_CK1
        use pm_kind, only: CKG => CK1
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#undef isCloseMean_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCloseMean_RK_ENABLED 1

#if RK5_ENABLED
    module procedure isCloseMean_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isCloseMean_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isCloseMean_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isCloseMean_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isCloseMean_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathCompare@routines.inc.F90"
    end procedure
#endif

#undef isCloseMean_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isCloseMean_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isClose_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines