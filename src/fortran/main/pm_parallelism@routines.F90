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
!>  This file contains procedure implementations of [pm_parallelism](@ref pm_parallelism).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_parallelism) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_val2str, only: getStr
    use pm_arraySpace, only: setLinSpace
    use pm_arrayResize, only: setResized
    use pm_distGeomCyclic, only: getGeomCyclicLogPMF

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setForkJoinScaling_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setForkJoinScaling_RK5
        use pm_kind, only: RKG => RK5
#include "pm_parallelism@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setForkJoinScaling_RK4
        use pm_kind, only: RKG => RK4
#include "pm_parallelism@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setForkJoinScaling_RK3
        use pm_kind, only: RKG => RK3
#include "pm_parallelism@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setForkJoinScaling_RK2
        use pm_kind, only: RKG => RK2
#include "pm_parallelism@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setForkJoinScaling_RK1
        use pm_kind, only: RKG => RK1
#include "pm_parallelism@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setForkJoinScaling_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines