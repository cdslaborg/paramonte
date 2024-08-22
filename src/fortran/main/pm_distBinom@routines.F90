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
!>  This file contains procedure implementations of [pm_distBinom](@ref pm_distBinom).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distBinom) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_mathBeta, only: setBetaInc
    use pm_mathBeta, only: getLogBeta

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getBinomLogPMF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getBinomLogPMF_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getBinomLogPMF_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getBinomLogPMF_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getBinomLogPMF_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getBinomLogPMF_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getBinomLogPMF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setBinomLogPMF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBinomLogPMF_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBinomLogPMF_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBinomLogPMF_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBinomLogPMF_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBinomLogPMF_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setBinomLogPMF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getBinomCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getBinomCDF_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getBinomCDF_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getBinomCDF_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getBinomCDF_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getBinomCDF_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getBinomCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setBinomCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBinomCDF_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBinomCDF_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBinomCDF_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBinomCDF_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBinomCDF_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distBinom@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setBinomCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines