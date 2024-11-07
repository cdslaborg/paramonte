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
!>  This file contains procedure implementations of [pm_knn](@ref pm_knn).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Wednesday 5:43 PM, December 25, 2013, Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_knn) routines

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define DO_CONCURRENT(I,IBEG,IEND) do I = IBEG, IEND
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define DO_CONCURRENT(I,IBEG,IEND) do concurrent(I = IBEG : IEND)
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_arraySort, only: setSorted
    use pm_arraySelect, only: setSelected
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setKnnSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Val_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setKnnSortedVal_RK5
        use pm_kind, only: TKG => RK5
#include "pm_knn@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setKnnSortedVal_RK4
        use pm_kind, only: TKG => RK4
#include "pm_knn@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setKnnSortedVal_RK3
        use pm_kind, only: TKG => RK3
#include "pm_knn@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setKnnSortedVal_RK2
        use pm_kind, only: TKG => RK2
#include "pm_knn@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setKnnSortedVal_RK1
        use pm_kind, only: TKG => RK1
#include "pm_knn@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Val_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Kth_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setKnnSortedKth_RK5
        use pm_kind, only: TKG => RK5
#include "pm_knn@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setKnnSortedKth_RK4
        use pm_kind, only: TKG => RK4
#include "pm_knn@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setKnnSortedKth_RK3
        use pm_kind, only: TKG => RK3
#include "pm_knn@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setKnnSortedKth_RK2
        use pm_kind, only: TKG => RK2
#include "pm_knn@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setKnnSortedKth_RK1
        use pm_kind, only: TKG => RK1
#include "pm_knn@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Kth_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ind_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setKnnSortedInd_RK5
        use pm_kind, only: TKG => RK5
#include "pm_knn@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setKnnSortedInd_RK4
        use pm_kind, only: TKG => RK4
#include "pm_knn@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setKnnSortedInd_RK3
        use pm_kind, only: TKG => RK3
#include "pm_knn@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setKnnSortedInd_RK2
        use pm_kind, only: TKG => RK2
#include "pm_knn@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setKnnSortedInd_RK1
        use pm_kind, only: TKG => RK1
#include "pm_knn@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ind_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setKnnSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION
#undef DO_CONCURRENT

end submodule routines