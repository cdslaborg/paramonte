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
!>  This file contains procedure implementations of [pm_distUnifPar](@ref pm_distUnifPar).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distUnifPar) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_matrixDet, only: setMatDetSqrtLog, uppDia, transHerm
    use pm_distUnif, only: setUnifRand
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getUnifParLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Cub_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifCubLogPDF_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifCubLogPDF_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifCubLogPDF_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifCubLogPDF_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifCubLogPDF_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Cub_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Rec_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRecLogPDF_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRecLogPDF_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRecLogPDF_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRecLogPDF_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRecLogPDF_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Rec_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Par_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifParLogPDF_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifParLogPDF_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifParLogPDF_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifParLogPDF_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifParLogPDF_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Par_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getUnifParLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getUnifParRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Cub_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifCubRandDU_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifCubRandDU_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifCubRandDU_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifCubRandDU_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifCubRandDU_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifCubRandLU_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifCubRandLU_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifCubRandLU_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifCubRandLU_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifCubRandLU_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Cub_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Rec_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRecRandDU_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRecRandDU_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRecRandDU_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRecRandDU_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRecRandDU_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifRecRandLU_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifRecRandLU_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifRecRandLU_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifRecRandLU_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifRecRandLU_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifParRandDU_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifParRandDU_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifParRandDU_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifParRandDU_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifParRandDU_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifParRandLU_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifParRandLU_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifParRandLU_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifParRandLU_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifParRandLU_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Rec_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getUnifParRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setUnifParRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Cub_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifCubRandDU_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifCubRandDU_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifCubRandDU_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifCubRandDU_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifCubRandDU_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifCubRandLU_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifCubRandLU_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifCubRandLU_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifCubRandLU_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifCubRandLU_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Cub_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Rec_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRecRandDU_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRecRandDU_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRecRandDU_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRecRandDU_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRecRandDU_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifRecRandLU_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifRecRandLU_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifRecRandLU_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifRecRandLU_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifRecRandLU_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Rec_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Par_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifParRandDU_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifParRandDU_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifParRandDU_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifParRandDU_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifParRandDU_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUnifParRandLU_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUnifParRandLU_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUnifParRandLU_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUnifParRandLU_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUnifParRandLU_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distUnifPar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Par_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setUnifParRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines
