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
!>  This file contains procedure implementations of [test_pm_arrayResize](@ref test_pm_arrayResize).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_arrayResize) routines ! LCOV_EXCL_LINE

    !   \bug
    !   There is a viscous Intel compiler 2022 bug where the appearance of the following `use` statements 
    !   in the body of the implementation include file `test_pm_arrayResize@routines.inc.F90` leads to various
    !   mistakes in parsing and preprocessing the contents of the include file.<br>
    !   The threshold for the maximum number of `use` statements within the entire submodule appears to be 
    !   about `55`, because activating more than 55 procedures of the submodule 
    !   leads to compilation failures due syntax parsing mistakes by the Intel compiler.<br>
    use pm_distUnif, only: getUnifRand
    use pm_distUnif, only: setUnifRand
    use pm_arrayInit, only: setCoreHalo
    use pm_io, only: display_type
    use pm_val2str, only: getStr

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D0_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setResized_D0_SK5_1
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setResized_D0_SK4_1
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setResized_D0_SK3_1
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setResized_D0_SK2_1
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setResized_D0_SK1_1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef setResized_D0_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setResized_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D1_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setResized_D1_SK5_1
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setResized_D1_SK4_1
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setResized_D1_SK3_1
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setResized_D1_SK2_1
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setResized_D1_SK1_1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef setResized_D1_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D1_IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setResized_D1_IK5_1
        use pm_kind, only: IKC => IK5
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setResized_D1_IK4_1
        use pm_kind, only: IKC => IK4
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setResized_D1_IK3_1
        use pm_kind, only: IKC => IK3
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setResized_D1_IK2_1
        use pm_kind, only: IKC => IK2
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setResized_D1_IK1_1
        use pm_kind, only: IKC => IK1
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef setResized_D1_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D1_LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setResized_D1_LK5_1
        use pm_kind, only: LKC => LK5
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setResized_D1_LK4_1
        use pm_kind, only: LKC => LK4
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setResized_D1_LK3_1
        use pm_kind, only: LKC => LK3
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setResized_D1_LK2_1
        use pm_kind, only: LKC => LK2
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setResized_D1_LK1_1
        use pm_kind, only: LKC => LK1
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef setResized_D1_LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D1_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setResized_D1_CK5_1
        use pm_kind, only: CKC => CK5
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setResized_D1_CK4_1
        use pm_kind, only: CKC => CK4
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setResized_D1_CK3_1
        use pm_kind, only: CKC => CK3
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setResized_D1_CK2_1
        use pm_kind, only: CKC => CK2
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setResized_D1_CK1_1
        use pm_kind, only: CKC => CK1
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef setResized_D1_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D1_RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setResized_D1_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setResized_D1_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setResized_D1_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setResized_D1_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setResized_D1_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef setResized_D1_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setResized_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D2_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setResized_D2_SK5_1
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setResized_D2_SK4_1
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setResized_D2_SK3_1
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setResized_D2_SK2_1
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setResized_D2_SK1_1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef setResized_D2_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D2_IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setResized_D2_IK5_1
        use pm_kind, only: IKC => IK5
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setResized_D2_IK4_1
        use pm_kind, only: IKC => IK4
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setResized_D2_IK3_1
        use pm_kind, only: IKC => IK3
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setResized_D2_IK2_1
        use pm_kind, only: IKC => IK2
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setResized_D2_IK1_1
        use pm_kind, only: IKC => IK1
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef setResized_D2_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D2_LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setResized_D2_LK5_1
        use pm_kind, only: LKC => LK5
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setResized_D2_LK4_1
        use pm_kind, only: LKC => LK4
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setResized_D2_LK3_1
        use pm_kind, only: LKC => LK3
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setResized_D2_LK2_1
        use pm_kind, only: LKC => LK2
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setResized_D2_LK1_1
        use pm_kind, only: LKC => LK1
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef setResized_D2_LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D2_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setResized_D2_CK5_1
        use pm_kind, only: CKC => CK5
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setResized_D2_CK4_1
        use pm_kind, only: CKC => CK4
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setResized_D2_CK3_1
        use pm_kind, only: CKC => CK3
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setResized_D2_CK2_1
        use pm_kind, only: CKC => CK2
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setResized_D2_CK1_1
        use pm_kind, only: CKC => CK1
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef setResized_D2_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D2_RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setResized_D2_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setResized_D2_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setResized_D2_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setResized_D2_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setResized_D2_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef setResized_D2_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setResized_D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D3_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setResized_D3_SK5_1
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setResized_D3_SK4_1
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setResized_D3_SK3_1
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setResized_D3_SK2_1
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setResized_D3_SK1_1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef setResized_D3_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D3_IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setResized_D3_IK5_1
        use pm_kind, only: IKC => IK5
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setResized_D3_IK4_1
        use pm_kind, only: IKC => IK4
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setResized_D3_IK3_1
        use pm_kind, only: IKC => IK3
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setResized_D3_IK2_1
        use pm_kind, only: IKC => IK2
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setResized_D3_IK1_1
        use pm_kind, only: IKC => IK1
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef setResized_D3_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D3_LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setResized_D3_LK5_1
        use pm_kind, only: LKC => LK5
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setResized_D3_LK4_1
        use pm_kind, only: LKC => LK4
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setResized_D3_LK3_1
        use pm_kind, only: LKC => LK3
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setResized_D3_LK2_1
        use pm_kind, only: LKC => LK2
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setResized_D3_LK1_1
        use pm_kind, only: LKC => LK1
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef setResized_D3_LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D3_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setResized_D3_CK5_1
        use pm_kind, only: CKC => CK5
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setResized_D3_CK4_1
        use pm_kind, only: CKC => CK4
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setResized_D3_CK3_1
        use pm_kind, only: CKC => CK3
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setResized_D3_CK2_1
        use pm_kind, only: CKC => CK2
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setResized_D3_CK1_1
        use pm_kind, only: CKC => CK1
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef setResized_D3_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setResized_D3_RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setResized_D3_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setResized_D3_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setResized_D3_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setResized_D3_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setResized_D3_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef setResized_D3_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setResized_D3_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setResized_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines