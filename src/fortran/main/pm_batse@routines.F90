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
!>  This file contains procedure implementations of [pm_batse](@ref pm_batse).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_batse) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_err, only: setAsserted
    use pm_arrayUnique, only: getUnique
    use pm_arraySort, only: isAscending
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_quadPack, only: getQuadErr, GK21, weps
    use pm_distGamma, only: setGammaCDF

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCorrectionLogEffPPF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCorrectionLogEffPPF_RK5
        use pm_kind, only: RKC => RK5
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCorrectionLogEffPPF_RK4
        use pm_kind, only: RKC => RK4
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCorrectionLogEffPPF_RK3
        use pm_kind, only: RKC => RK3
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCorrectionLogEffPPF_RK2
        use pm_kind, only: RKC => RK2
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCorrectionLogEffPPF_RK1
        use pm_kind, only: RKC => RK1
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCorrectionLogEffPPF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogEffPPF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogEffPPF_RK5
        use pm_kind, only: RKC => RK5
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogEffPPF_RK4
        use pm_kind, only: RKC => RK4
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogEffPPF_RK3
        use pm_kind, only: RKC => RK3
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogEffPPF_RK2
        use pm_kind, only: RKC => RK2
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogEffPPF_RK1
        use pm_kind, only: RKC => RK1
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogEffPPF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogPbol_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogPbol_RK5
        use pm_kind, only: RKC => RK5
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogPbol_RK4
        use pm_kind, only: RKC => RK4
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogPbol_RK3
        use pm_kind, only: RKC => RK3
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogPbol_RK2
        use pm_kind, only: RKC => RK2
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogPbol_RK1
        use pm_kind, only: RKC => RK1
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogPbol_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogPF53_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogPF53_RK5
        use pm_kind, only: RKC => RK5
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogPF53_RK4
        use pm_kind, only: RKC => RK4
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogPF53_RK3
        use pm_kind, only: RKC => RK3
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogPF53_RK2
        use pm_kind, only: RKC => RK2
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogPF53_RK1
        use pm_kind, only: RKC => RK1
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogPF53_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLog10PF53_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLog10PF53_RK5
        use pm_kind, only: RKC => RK5
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLog10PF53_RK4
        use pm_kind, only: RKC => RK4
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLog10PF53_RK3
        use pm_kind, only: RKC => RK3
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLog10PF53_RK2
        use pm_kind, only: RKC => RK2
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLog10PF53_RK1
        use pm_kind, only: RKC => RK1
#include "pm_batse@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLog10PF53_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines ! LCOV_EXCL_LINE