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
!>  This file contains procedure implementations of [pm_distChol](@ref pm_distChol).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distChol) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_distUnif, only: setUnifRand
    use pm_matrixInit, only: setMatInit, upp, low

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCholRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getCholRandRNGD_CK5
        use pm_kind, only: TKG => CK5
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCholRandRNGD_CK4
        use pm_kind, only: TKG => CK4
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCholRandRNGD_CK3
        use pm_kind, only: TKG => CK3
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCholRandRNGD_CK2
        use pm_kind, only: TKG => CK2
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCholRandRNGD_CK1
        use pm_kind, only: TKG => CK1
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCholRandRNGD_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCholRandRNGD_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCholRandRNGD_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCholRandRNGD_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCholRandRNGD_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCholRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCholRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCholRandRNGF_CK5
        use pm_kind, only: TKG => CK5
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCholRandRNGF_CK4
        use pm_kind, only: TKG => CK4
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCholRandRNGF_CK3
        use pm_kind, only: TKG => CK3
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCholRandRNGF_CK2
        use pm_kind, only: TKG => CK2
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCholRandRNGF_CK1
        use pm_kind, only: TKG => CK1
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCholRandRNGF_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCholRandRNGF_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCholRandRNGF_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCholRandRNGF_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCholRandRNGF_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distChol@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCholRandRNGX_CK5
        use pm_kind, only: TKG => CK5
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCholRandRNGX_CK4
        use pm_kind, only: TKG => CK4
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCholRandRNGX_CK3
        use pm_kind, only: TKG => CK3
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCholRandRNGX_CK2
        use pm_kind, only: TKG => CK2
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCholRandRNGX_CK1
        use pm_kind, only: TKG => CK1
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCholRandRNGX_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCholRandRNGX_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCholRandRNGX_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCholRandRNGX_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCholRandRNGX_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distChol@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCholRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines