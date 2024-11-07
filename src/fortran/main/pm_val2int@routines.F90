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
!>  This file contains procedure implementations of [val2int_pmod](@ref pm_val2int).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_val2int) routines ! LCOV_EXCL_LINE

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getInt_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getIntDef_LK5
        use pm_kind, only: IKG => IK, LKG => LK5
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getIntDef_LK4
        use pm_kind, only: IKG => IK, LKG => LK4
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getIntDef_LK3
        use pm_kind, only: IKG => IK, LKG => LK3
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getIntDef_LK2
        use pm_kind, only: IKG => IK, LKG => LK2
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getIntDef_LK1
        use pm_kind, only: IKG => IK, LKG => LK1
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getIntDef_SK5
        use pm_kind, only: IKG => IK, SKG => SK5
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getIntDef_SK4
        use pm_kind, only: IKG => IK, SKG => SK4
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getIntDef_SK3
        use pm_kind, only: IKG => IK, SKG => SK3
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getIntDef_SK2
        use pm_kind, only: IKG => IK, SKG => SK2
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getIntDef_SK1
        use pm_kind, only: IKG => IK, SKG => SK1
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getInt_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setInt_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED && LK5_ENABLED
    module procedure setIntDef_IK5_LK5
        use pm_kind, only: IKG => IK5, LKG => LK5
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && LK4_ENABLED
    module procedure setIntDef_IK5_LK4
        use pm_kind, only: IKG => IK5, LKG => LK4
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && LK3_ENABLED
    module procedure setIntDef_IK5_LK3
        use pm_kind, only: IKG => IK5, LKG => LK3
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && LK2_ENABLED
    module procedure setIntDef_IK5_LK2
        use pm_kind, only: IKG => IK5, LKG => LK2
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && LK1_ENABLED
    module procedure setIntDef_IK5_LK1
        use pm_kind, only: IKG => IK5, LKG => LK1
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK4_ENABLED && LK5_ENABLED
    module procedure setIntDef_IK4_LK5
        use pm_kind, only: IKG => IK4, LKG => LK5
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && LK4_ENABLED
    module procedure setIntDef_IK4_LK4
        use pm_kind, only: IKG => IK4, LKG => LK4
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && LK3_ENABLED
    module procedure setIntDef_IK4_LK3
        use pm_kind, only: IKG => IK4, LKG => LK3
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && LK2_ENABLED
    module procedure setIntDef_IK4_LK2
        use pm_kind, only: IKG => IK4, LKG => LK2
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && LK1_ENABLED
    module procedure setIntDef_IK4_LK1
        use pm_kind, only: IKG => IK4, LKG => LK1
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK3_ENABLED && LK5_ENABLED
    module procedure setIntDef_IK3_LK5
        use pm_kind, only: IKG => IK3, LKG => LK5
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && LK4_ENABLED
    module procedure setIntDef_IK3_LK4
        use pm_kind, only: IKG => IK3, LKG => LK4
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && LK3_ENABLED
    module procedure setIntDef_IK3_LK3
        use pm_kind, only: IKG => IK3, LKG => LK3
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && LK2_ENABLED
    module procedure setIntDef_IK3_LK2
        use pm_kind, only: IKG => IK3, LKG => LK2
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && LK1_ENABLED
    module procedure setIntDef_IK3_LK1
        use pm_kind, only: IKG => IK3, LKG => LK1
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK2_ENABLED && LK5_ENABLED
    module procedure setIntDef_IK2_LK5
        use pm_kind, only: IKG => IK2, LKG => LK5
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && LK4_ENABLED
    module procedure setIntDef_IK2_LK4
        use pm_kind, only: IKG => IK2, LKG => LK4
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && LK3_ENABLED
    module procedure setIntDef_IK2_LK3
        use pm_kind, only: IKG => IK2, LKG => LK3
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && LK2_ENABLED
    module procedure setIntDef_IK2_LK2
        use pm_kind, only: IKG => IK2, LKG => LK2
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && LK1_ENABLED
    module procedure setIntDef_IK2_LK1
        use pm_kind, only: IKG => IK2, LKG => LK1
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK1_ENABLED && LK5_ENABLED
    module procedure setIntDef_IK1_LK5
        use pm_kind, only: IKG => IK1, LKG => LK5
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && LK4_ENABLED
    module procedure setIntDef_IK1_LK4
        use pm_kind, only: IKG => IK1, LKG => LK4
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && LK3_ENABLED
    module procedure setIntDef_IK1_LK3
        use pm_kind, only: IKG => IK1, LKG => LK3
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && LK2_ENABLED
    module procedure setIntDef_IK1_LK2
        use pm_kind, only: IKG => IK1, LKG => LK2
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && LK1_ENABLED
    module procedure setIntDef_IK1_LK1
        use pm_kind, only: IKG => IK1, LKG => LK1
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED && SK5_ENABLED
    module procedure setIntDef_IK5_SK5
        use pm_kind, only: IKG => IK5, SKG => SK5
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && SK4_ENABLED
    module procedure setIntDef_IK5_SK4
        use pm_kind, only: IKG => IK5, SKG => SK4
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && SK3_ENABLED
    module procedure setIntDef_IK5_SK3
        use pm_kind, only: IKG => IK5, SKG => SK3
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && SK2_ENABLED
    module procedure setIntDef_IK5_SK2
        use pm_kind, only: IKG => IK5, SKG => SK2
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && SK1_ENABLED
    module procedure setIntDef_IK5_SK1
        use pm_kind, only: IKG => IK5, SKG => SK1
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK4_ENABLED && SK5_ENABLED
    module procedure setIntDef_IK4_SK5
        use pm_kind, only: IKG => IK4, SKG => SK5
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && SK4_ENABLED
    module procedure setIntDef_IK4_SK4
        use pm_kind, only: IKG => IK4, SKG => SK4
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && SK3_ENABLED
    module procedure setIntDef_IK4_SK3
        use pm_kind, only: IKG => IK4, SKG => SK3
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && SK2_ENABLED
    module procedure setIntDef_IK4_SK2
        use pm_kind, only: IKG => IK4, SKG => SK2
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && SK1_ENABLED
    module procedure setIntDef_IK4_SK1
        use pm_kind, only: IKG => IK4, SKG => SK1
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK3_ENABLED && SK5_ENABLED
    module procedure setIntDef_IK3_SK5
        use pm_kind, only: IKG => IK3, SKG => SK5
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && SK4_ENABLED
    module procedure setIntDef_IK3_SK4
        use pm_kind, only: IKG => IK3, SKG => SK4
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && SK3_ENABLED
    module procedure setIntDef_IK3_SK3
        use pm_kind, only: IKG => IK3, SKG => SK3
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && SK2_ENABLED
    module procedure setIntDef_IK3_SK2
        use pm_kind, only: IKG => IK3, SKG => SK2
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && SK1_ENABLED
    module procedure setIntDef_IK3_SK1
        use pm_kind, only: IKG => IK3, SKG => SK1
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK2_ENABLED && SK5_ENABLED
    module procedure setIntDef_IK2_SK5
        use pm_kind, only: IKG => IK2, SKG => SK5
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && SK4_ENABLED
    module procedure setIntDef_IK2_SK4
        use pm_kind, only: IKG => IK2, SKG => SK4
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && SK3_ENABLED
    module procedure setIntDef_IK2_SK3
        use pm_kind, only: IKG => IK2, SKG => SK3
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && SK2_ENABLED
    module procedure setIntDef_IK2_SK2
        use pm_kind, only: IKG => IK2, SKG => SK2
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && SK1_ENABLED
    module procedure setIntDef_IK2_SK1
        use pm_kind, only: IKG => IK2, SKG => SK1
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK1_ENABLED && SK5_ENABLED
    module procedure setIntDef_IK1_SK5
        use pm_kind, only: IKG => IK1, SKG => SK5
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && SK4_ENABLED
    module procedure setIntDef_IK1_SK4
        use pm_kind, only: IKG => IK1, SKG => SK4
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && SK3_ENABLED
    module procedure setIntDef_IK1_SK3
        use pm_kind, only: IKG => IK1, SKG => SK3
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && SK2_ENABLED
    module procedure setIntDef_IK1_SK2
        use pm_kind, only: IKG => IK1, SKG => SK2
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && SK1_ENABLED
    module procedure setIntDef_IK1_SK1
        use pm_kind, only: IKG => IK1, SKG => SK1
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Err_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED && SK5_ENABLED
    module procedure setIntErr_IK5_SK5
        use pm_kind, only: IKG => IK5, SKG => SK5
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && SK4_ENABLED
    module procedure setIntErr_IK5_SK4
        use pm_kind, only: IKG => IK5, SKG => SK4
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && SK3_ENABLED
    module procedure setIntErr_IK5_SK3
        use pm_kind, only: IKG => IK5, SKG => SK3
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && SK2_ENABLED
    module procedure setIntErr_IK5_SK2
        use pm_kind, only: IKG => IK5, SKG => SK2
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && SK1_ENABLED
    module procedure setIntErr_IK5_SK1
        use pm_kind, only: IKG => IK5, SKG => SK1
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK4_ENABLED && SK5_ENABLED
    module procedure setIntErr_IK4_SK5
        use pm_kind, only: IKG => IK4, SKG => SK5
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && SK4_ENABLED
    module procedure setIntErr_IK4_SK4
        use pm_kind, only: IKG => IK4, SKG => SK4
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && SK3_ENABLED
    module procedure setIntErr_IK4_SK3
        use pm_kind, only: IKG => IK4, SKG => SK3
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && SK2_ENABLED
    module procedure setIntErr_IK4_SK2
        use pm_kind, only: IKG => IK4, SKG => SK2
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && SK1_ENABLED
    module procedure setIntErr_IK4_SK1
        use pm_kind, only: IKG => IK4, SKG => SK1
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK3_ENABLED && SK5_ENABLED
    module procedure setIntErr_IK3_SK5
        use pm_kind, only: IKG => IK3, SKG => SK5
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && SK4_ENABLED
    module procedure setIntErr_IK3_SK4
        use pm_kind, only: IKG => IK3, SKG => SK4
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && SK3_ENABLED
    module procedure setIntErr_IK3_SK3
        use pm_kind, only: IKG => IK3, SKG => SK3
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && SK2_ENABLED
    module procedure setIntErr_IK3_SK2
        use pm_kind, only: IKG => IK3, SKG => SK2
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && SK1_ENABLED
    module procedure setIntErr_IK3_SK1
        use pm_kind, only: IKG => IK3, SKG => SK1
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK2_ENABLED && SK5_ENABLED
    module procedure setIntErr_IK2_SK5
        use pm_kind, only: IKG => IK2, SKG => SK5
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && SK4_ENABLED
    module procedure setIntErr_IK2_SK4
        use pm_kind, only: IKG => IK2, SKG => SK4
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && SK3_ENABLED
    module procedure setIntErr_IK2_SK3
        use pm_kind, only: IKG => IK2, SKG => SK3
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && SK2_ENABLED
    module procedure setIntErr_IK2_SK2
        use pm_kind, only: IKG => IK2, SKG => SK2
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && SK1_ENABLED
    module procedure setIntErr_IK2_SK1
        use pm_kind, only: IKG => IK2, SKG => SK1
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK1_ENABLED && SK5_ENABLED
    module procedure setIntErr_IK1_SK5
        use pm_kind, only: IKG => IK1, SKG => SK5
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && SK4_ENABLED
    module procedure setIntErr_IK1_SK4
        use pm_kind, only: IKG => IK1, SKG => SK4
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && SK3_ENABLED
    module procedure setIntErr_IK1_SK3
        use pm_kind, only: IKG => IK1, SKG => SK3
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && SK2_ENABLED
    module procedure setIntErr_IK1_SK2
        use pm_kind, only: IKG => IK1, SKG => SK2
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && SK1_ENABLED
    module procedure setIntErr_IK1_SK1
        use pm_kind, only: IKG => IK1, SKG => SK1
#include "pm_val2int@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Err_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setInt_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
