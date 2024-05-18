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
!>  This file contains procedure implementations of [pm_arrayChange](@ref pm_arrayChange).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_arrayChange) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_arrayRange, only: getRange
    use pm_arrayChoice, only: setChoice
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getChange_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unif_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getChangeUnifRNGD_SK5
        use pm_kind, only: SKC => SK5
        use pm_distUnif, only: rng => rngf
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getChangeUnifRNGD_SK4
        use pm_kind, only: SKC => SK4
        use pm_distUnif, only: rng => rngf
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getChangeUnifRNGD_SK3
        use pm_kind, only: SKC => SK3
        use pm_distUnif, only: rng => rngf
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getChangeUnifRNGD_SK2
        use pm_kind, only: SKC => SK2
        use pm_distUnif, only: rng => rngf
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getChangeUnifRNGD_SK1
        use pm_kind, only: SKC => SK1
        use pm_distUnif, only: rng => rngf
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getChangeUnifRNGD_IK5
        use pm_kind, only: IKC => IK5
        use pm_distUnif, only: rng => rngf
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getChangeUnifRNGD_IK4
        use pm_kind, only: IKC => IK4
        use pm_distUnif, only: rng => rngf
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getChangeUnifRNGD_IK3
        use pm_kind, only: IKC => IK3
        use pm_distUnif, only: rng => rngf
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getChangeUnifRNGD_IK2
        use pm_kind, only: IKC => IK2
        use pm_distUnif, only: rng => rngf
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getChangeUnifRNGD_IK1
        use pm_kind, only: IKC => IK1
        use pm_distUnif, only: rng => rngf
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getChangeUnifRNGD_RK5
        use pm_kind, only: RKC => RK5
        use pm_distUnif, only: rng => rngf
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getChangeUnifRNGD_RK4
        use pm_kind, only: RKC => RK4
        use pm_distUnif, only: rng => rngf
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getChangeUnifRNGD_RK3
        use pm_kind, only: RKC => RK3
        use pm_distUnif, only: rng => rngf
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getChangeUnifRNGD_RK2
        use pm_kind, only: RKC => RK2
        use pm_distUnif, only: rng => rngf
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getChangeUnifRNGD_RK1
        use pm_kind, only: RKC => RK1
        use pm_distUnif, only: rng => rngf
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Unif_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getChange_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setChange_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unif_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setChangeUnifRNGF_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setChangeUnifRNGF_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setChangeUnifRNGF_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setChangeUnifRNGF_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setChangeUnifRNGF_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setChangeUnifRNGF_IK5
        use pm_kind, only: IKC => IK5
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setChangeUnifRNGF_IK4
        use pm_kind, only: IKC => IK4
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setChangeUnifRNGF_IK3
        use pm_kind, only: IKC => IK3
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setChangeUnifRNGF_IK2
        use pm_kind, only: IKC => IK2
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setChangeUnifRNGF_IK1
        use pm_kind, only: IKC => IK1
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setChangeUnifRNGF_RK5
        use pm_kind, only: RKC => RK5
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setChangeUnifRNGF_RK4
        use pm_kind, only: RKC => RK4
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setChangeUnifRNGF_RK3
        use pm_kind, only: RKC => RK3
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setChangeUnifRNGF_RK2
        use pm_kind, only: RKC => RK2
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setChangeUnifRNGF_RK1
        use pm_kind, only: RKC => RK1
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setChangeUnifRNGX_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setChangeUnifRNGX_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setChangeUnifRNGX_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setChangeUnifRNGX_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setChangeUnifRNGX_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setChangeUnifRNGX_IK5
        use pm_kind, only: IKC => IK5
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setChangeUnifRNGX_IK4
        use pm_kind, only: IKC => IK4
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setChangeUnifRNGX_IK3
        use pm_kind, only: IKC => IK3
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setChangeUnifRNGX_IK2
        use pm_kind, only: IKC => IK2
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setChangeUnifRNGX_IK1
        use pm_kind, only: IKC => IK1
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setChangeUnifRNGX_RK5
        use pm_kind, only: RKC => RK5
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setChangeUnifRNGX_RK4
        use pm_kind, only: RKC => RK4
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setChangeUnifRNGX_RK3
        use pm_kind, only: RKC => RK3
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setChangeUnifRNGX_RK2
        use pm_kind, only: RKC => RK2
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setChangeUnifRNGX_RK1
        use pm_kind, only: RKC => RK1
#include "pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Unif_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setChange_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines
