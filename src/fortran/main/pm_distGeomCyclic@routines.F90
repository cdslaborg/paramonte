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
!>  This file contains procedure implementations of [pm_distGeomCyclic](@ref pm_distGeomCyclic).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distGeomCyclic) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_math1mexp, only: get1mexp
    use pm_distGeom, only: setGeomRand
    use pm_optimization, only: getMinBrent
    use pm_optimization, only: setMinBrent
    use pm_optimization, only: isFailedMinPowell

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGeomCyclicLogPMF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGeomCyclicLogPMF_D0_IK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGeomCyclicLogPMF_D0_IK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGeomCyclicLogPMF_D0_IK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGeomCyclicLogPMF_D0_IK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGeomCyclicLogPMF_D0_IK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGeomCyclicLogPMF_D0_RK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGeomCyclicLogPMF_D0_RK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGeomCyclicLogPMF_D0_RK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGeomCyclicLogPMF_D0_RK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGeomCyclicLogPMF_D0_RK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGeomCyclicLogPMF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGeomCyclicLogPMF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicLogPMFDef_D0_IK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicLogPMFDef_D0_IK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicLogPMFDef_D0_IK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicLogPMFDef_D0_IK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicLogPMFDef_D0_IK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicLogPMFDef_D0_RK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicLogPMFDef_D0_RK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicLogPMFDef_D0_RK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicLogPMFDef_D0_RK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicLogPMFDef_D0_RK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Log_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicLogPMFLog_D0_IK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicLogPMFLog_D0_IK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicLogPMFLog_D0_IK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicLogPMFLog_D0_IK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicLogPMFLog_D0_IK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicLogPMFLog_D0_RK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicLogPMFLog_D0_RK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicLogPMFLog_D0_RK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicLogPMFLog_D0_RK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicLogPMFLog_D0_RK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Log_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicLogPMFDef_D1_IK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicLogPMFDef_D1_IK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicLogPMFDef_D1_IK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicLogPMFDef_D1_IK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicLogPMFDef_D1_IK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicLogPMFDef_D1_RK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicLogPMFDef_D1_RK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicLogPMFDef_D1_RK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicLogPMFDef_D1_RK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicLogPMFDef_D1_RK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Log_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicLogPMFLog_D1_IK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicLogPMFLog_D1_IK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicLogPMFLog_D1_IK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicLogPMFLog_D1_IK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicLogPMFLog_D1_IK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicLogPMFLog_D1_RK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicLogPMFLog_D1_RK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicLogPMFLog_D1_RK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicLogPMFLog_D1_RK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicLogPMFLog_D1_RK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Log_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGeomCyclicLogPMF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGeomCyclicLogCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGeomCyclicLogCDF_D0_IK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGeomCyclicLogCDF_D0_IK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGeomCyclicLogCDF_D0_IK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGeomCyclicLogCDF_D0_IK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGeomCyclicLogCDF_D0_IK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGeomCyclicLogCDF_D0_RK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGeomCyclicLogCDF_D0_RK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGeomCyclicLogCDF_D0_RK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGeomCyclicLogCDF_D0_RK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGeomCyclicLogCDF_D0_RK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGeomCyclicLogCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGeomCyclicLogCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicLogCDFDef_D0_IK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicLogCDFDef_D0_IK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicLogCDFDef_D0_IK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicLogCDFDef_D0_IK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicLogCDFDef_D0_IK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicLogCDFDef_D0_RK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicLogCDFDef_D0_RK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicLogCDFDef_D0_RK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicLogCDFDef_D0_RK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicLogCDFDef_D0_RK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Log_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicLogCDFLog_D0_IK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicLogCDFLog_D0_IK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicLogCDFLog_D0_IK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicLogCDFLog_D0_IK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicLogCDFLog_D0_IK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicLogCDFLog_D0_RK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicLogCDFLog_D0_RK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicLogCDFLog_D0_RK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicLogCDFLog_D0_RK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicLogCDFLog_D0_RK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Log_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicLogCDFDef_D1_IK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicLogCDFDef_D1_IK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicLogCDFDef_D1_IK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicLogCDFDef_D1_IK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicLogCDFDef_D1_IK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicLogCDFDef_D1_RK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicLogCDFDef_D1_RK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicLogCDFDef_D1_RK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicLogCDFDef_D1_RK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicLogCDFDef_D1_RK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Log_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicLogCDFLog_D1_IK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicLogCDFLog_D1_IK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicLogCDFLog_D1_IK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicLogCDFLog_D1_IK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicLogCDFLog_D1_IK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicLogCDFLog_D1_RK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicLogCDFLog_D1_RK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicLogCDFLog_D1_RK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicLogCDFLog_D1_RK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicLogCDFLog_D1_RK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Log_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGeomCyclicLogCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGeomCyclicRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGeomCyclicRand_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGeomCyclicRand_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGeomCyclicRand_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGeomCyclicRand_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGeomCyclicRand_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGeomCyclicRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGeomCyclicRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicRandRNGD_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicRandRNGD_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicRandRNGD_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicRandRNGD_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicRandRNGD_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicRandRNGF_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicRandRNGF_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicRandRNGF_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicRandRNGF_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicRandRNGF_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
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

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicRandRNGX_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicRandRNGX_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicRandRNGX_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicRandRNGX_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicRandRNGX_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicRandRNGD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicRandRNGD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicRandRNGD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicRandRNGD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicRandRNGD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicRandRNGF_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicRandRNGF_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicRandRNGF_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicRandRNGF_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicRandRNGF_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
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

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCyclicRandRNGX_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCyclicRandRNGX_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCyclicRandRNGX_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCyclicRandRNGX_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCyclicRandRNGX_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGeomCyclicRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isFailedGeomCyclicFit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isFailedGeomCyclicFit_IK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isFailedGeomCyclicFit_IK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isFailedGeomCyclicFit_IK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isFailedGeomCyclicFit_IK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isFailedGeomCyclicFit_IK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isFailedGeomCyclicFit_RK_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isFailedGeomCyclicFit_RK_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isFailedGeomCyclicFit_RK_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isFailedGeomCyclicFit_RK_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isFailedGeomCyclicFit_RK_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distGeomCyclic@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isFailedGeomCyclicFit_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines