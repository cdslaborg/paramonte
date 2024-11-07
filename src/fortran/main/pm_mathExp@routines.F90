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
!>  This file contains procedure implementations of [pm_mathExp](@ref pm_mathExp).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_mathExp) routines ! LCOV_EXCL_LINE

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

#define isIntPow_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arb_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure isIntPowArb_IK5
        use pm_kind, only: IKG => IK5
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure isIntPowArb_IK4
        use pm_kind, only: IKG => IK4
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure isIntPowArb_IK3
        use pm_kind, only: IKG => IK3
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure isIntPowArb_IK2
        use pm_kind, only: IKG => IK2
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure isIntPowArb_IK1
        use pm_kind, only: IKG => IK1
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arb_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure isIntPowDef_IK5
        use pm_kind, only: IKG => IK5
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure isIntPowDef_IK4
        use pm_kind, only: IKG => IK4
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure isIntPowDef_IK3
        use pm_kind, only: IKG => IK3
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure isIntPowDef_IK2
        use pm_kind, only: IKG => IK2
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure isIntPowDef_IK1
        use pm_kind, only: IKG => IK1
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isIntPow_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getExpNext_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getExpNext_IK5
        use pm_kind, only: IKG => IK5
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getExpNext_IK4
        use pm_kind, only: IKG => IK4
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getExpNext_IK3
        use pm_kind, only: IKG => IK3
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getExpNext_IK2
        use pm_kind, only: IKG => IK2
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getExpNext_IK1
        use pm_kind, only: IKG => IK1
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getExpNext_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExpNext_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExpNext_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExpNext_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExpNext_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getExpNext_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getExpPrev_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getExpPrev_IK5
        use pm_kind, only: IKG => IK5
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getExpPrev_IK4
        use pm_kind, only: IKG => IK4
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getExpPrev_IK3
        use pm_kind, only: IKG => IK3
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getExpPrev_IK2
        use pm_kind, only: IKG => IK2
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getExpPrev_IK1
        use pm_kind, only: IKG => IK1
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getExpPrev_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExpPrev_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExpPrev_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExpPrev_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExpPrev_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getExpPrev_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines