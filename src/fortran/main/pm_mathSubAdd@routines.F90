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
!>  This file contains procedure implementations of [pm_mathSubAdd](@ref pm_mathSubAdd).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_mathSubAdd) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_kind, only: SK

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getSubAdd_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getSubAddUnary_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getSubAddUnary_IK_ENABLED 1

#if IK5_ENABLED
    module procedure getSubAddUnary_IK5
        use pm_kind, only: IKC => IK5
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getSubAddUnary_IK4
        use pm_kind, only: IKC => IK4
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getSubAddUnary_IK3
        use pm_kind, only: IKC => IK3
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getSubAddUnary_IK2
        use pm_kind, only: IKC => IK2
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getSubAddUnary_IK1
        use pm_kind, only: IKC => IK1
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#undef getSubAddUnary_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getSubAddUnary_RK_ENABLED 1

#if RK5_ENABLED
    module procedure getSubAddUnary_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getSubAddUnary_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getSubAddUnary_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getSubAddUnary_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getSubAddUnary_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#undef getSubAddUnary_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getSubAddUnary_CK_ENABLED 1

#if CK5_ENABLED
    module procedure getSubAddUnary_CK5
        use pm_kind, only: CKC => CK5
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getSubAddUnary_CK4
        use pm_kind, only: CKC => CK4
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getSubAddUnary_CK3
        use pm_kind, only: CKC => CK3
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getSubAddUnary_CK2
        use pm_kind, only: CKC => CK2
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getSubAddUnary_CK1
        use pm_kind, only: CKC => CK1
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#undef getSubAddUnary_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getSubAddUnary_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getSubAddBinary_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getSubAddBinary_IK_ENABLED 1

#if IK5_ENABLED
    module procedure getSubAddBinary_IK5
        use pm_kind, only: IKC => IK5
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getSubAddBinary_IK4
        use pm_kind, only: IKC => IK4
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getSubAddBinary_IK3
        use pm_kind, only: IKC => IK3
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getSubAddBinary_IK2
        use pm_kind, only: IKC => IK2
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getSubAddBinary_IK1
        use pm_kind, only: IKC => IK1
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#undef getSubAddBinary_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getSubAddBinary_RK_ENABLED 1

#if RK5_ENABLED
    module procedure getSubAddBinary_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getSubAddBinary_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getSubAddBinary_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getSubAddBinary_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getSubAddBinary_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#undef getSubAddBinary_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getSubAddBinary_CK_ENABLED 1

#if CK5_ENABLED
    module procedure getSubAddBinary_CK5
        use pm_kind, only: CKC => CK5
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getSubAddBinary_CK4
        use pm_kind, only: CKC => CK4
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getSubAddBinary_CK3
        use pm_kind, only: CKC => CK3
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getSubAddBinary_CK2
        use pm_kind, only: CKC => CK2
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getSubAddBinary_CK1
        use pm_kind, only: CKC => CK1
#include "pm_mathSubAdd@routines.inc.F90"
    end procedure
#endif

#undef getSubAddBinary_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getSubAddBinary_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getSubAdd_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
