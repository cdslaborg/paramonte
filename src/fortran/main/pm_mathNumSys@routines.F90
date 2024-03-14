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
!>  This file contains procedure implementations of [pm_mathNumSys](@ref pm_mathNumSys).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Wednesday 03:29 AM, September 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_mathNumSys) routines ! LCOV_EXCL_LINE

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

#define getDecimal_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDecimal_IK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED && IK5_ENABLED
    module procedure getDecimal_SK5_IK5
        use pm_kind, only: SKC => SK5, IKC => IK5
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && IK4_ENABLED
    module procedure getDecimal_SK5_IK4
        use pm_kind, only: SKC => SK5, IKC => IK4
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && IK3_ENABLED
    module procedure getDecimal_SK5_IK3
        use pm_kind, only: SKC => SK5, IKC => IK3
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && IK2_ENABLED
    module procedure getDecimal_SK5_IK2
        use pm_kind, only: SKC => SK5, IKC => IK2
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && IK1_ENABLED
    module procedure getDecimal_SK5_IK1
        use pm_kind, only: SKC => SK5, IKC => IK1
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK4_ENABLED && IK5_ENABLED
    module procedure getDecimal_SK4_IK5
        use pm_kind, only: SKC => SK4, IKC => IK5
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && IK4_ENABLED
    module procedure getDecimal_SK4_IK4
        use pm_kind, only: SKC => SK4, IKC => IK4
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && IK3_ENABLED
    module procedure getDecimal_SK4_IK3
        use pm_kind, only: SKC => SK4, IKC => IK3
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && IK2_ENABLED
    module procedure getDecimal_SK4_IK2
        use pm_kind, only: SKC => SK4, IKC => IK2
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && IK1_ENABLED
    module procedure getDecimal_SK4_IK1
        use pm_kind, only: SKC => SK4, IKC => IK1
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK3_ENABLED && IK5_ENABLED
    module procedure getDecimal_SK3_IK5
        use pm_kind, only: SKC => SK3, IKC => IK5
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && IK4_ENABLED
    module procedure getDecimal_SK3_IK4
        use pm_kind, only: SKC => SK3, IKC => IK4
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && IK3_ENABLED
    module procedure getDecimal_SK3_IK3
        use pm_kind, only: SKC => SK3, IKC => IK3
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && IK2_ENABLED
    module procedure getDecimal_SK3_IK2
        use pm_kind, only: SKC => SK3, IKC => IK2
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && IK1_ENABLED
    module procedure getDecimal_SK3_IK1
        use pm_kind, only: SKC => SK3, IKC => IK1
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK2_ENABLED && IK5_ENABLED
    module procedure getDecimal_SK2_IK5
        use pm_kind, only: SKC => SK2, IKC => IK5
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && IK4_ENABLED
    module procedure getDecimal_SK2_IK4
        use pm_kind, only: SKC => SK2, IKC => IK4
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && IK3_ENABLED
    module procedure getDecimal_SK2_IK3
        use pm_kind, only: SKC => SK2, IKC => IK3
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && IK2_ENABLED
    module procedure getDecimal_SK2_IK2
        use pm_kind, only: SKC => SK2, IKC => IK2
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && IK1_ENABLED
    module procedure getDecimal_SK2_IK1
        use pm_kind, only: SKC => SK2, IKC => IK1
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK1_ENABLED && IK5_ENABLED
    module procedure getDecimal_SK1_IK5
        use pm_kind, only: SKC => SK1, IKC => IK5
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && IK4_ENABLED
    module procedure getDecimal_SK1_IK4
        use pm_kind, only: SKC => SK1, IKC => IK4
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && IK3_ENABLED
    module procedure getDecimal_SK1_IK3
        use pm_kind, only: SKC => SK1, IKC => IK3
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && IK2_ENABLED
    module procedure getDecimal_SK1_IK2
        use pm_kind, only: SKC => SK1, IKC => IK2
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && IK1_ENABLED
    module procedure getDecimal_SK1_IK1
        use pm_kind, only: SKC => SK1, IKC => IK1
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDecimal_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDecimal_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNumeral_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNumeral_IK_ENABLED 1

#if IK5_ENABLED
    module procedure getNumeral_IK5
        use pm_kind, only: SKC => SK, IKC => IK5
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getNumeral_IK4
        use pm_kind, only: SKC => SK, IKC => IK4
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getNumeral_IK3
        use pm_kind, only: SKC => SK, IKC => IK3
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getNumeral_IK2
        use pm_kind, only: SKC => SK, IKC => IK2
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getNumeral_IK1
        use pm_kind, only: SKC => SK, IKC => IK1
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#undef getNumeral_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNumeral_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCountDigit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCountDigit_IK_ENABLED 1

#if IK5_ENABLED
    module procedure getCountDigit_IK5
        use pm_kind, only: IKC => IK5
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getCountDigit_IK4
        use pm_kind, only: IKC => IK4
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getCountDigit_IK3
        use pm_kind, only: IKC => IK3
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getCountDigit_IK2
        use pm_kind, only: IKC => IK2
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getCountDigit_IK1
        use pm_kind, only: IKC => IK1
#include "pm_mathNumSys@routines.inc.F90"
    end procedure
#endif

#undef getCountDigit_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCountDigit_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
