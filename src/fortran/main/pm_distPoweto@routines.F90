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
!>  This file contains procedure implementations of [pm_distPoweto](@ref pm_distPoweto).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distPoweto) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
    use pm_arraySort, only: isAscending
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_distNegExp, only: getNegExpRand
    use pm_mathLogSubExp, only: getLogSubExp
    use pm_distPower, only: getPowerLogPDFNF, setPowerLogPDF, getPowerLogCDFNF, setPowerLogCDF, setPowerLogQuan, setPowerLogRand
    use pm_distPareto, only: getParetoLogPDFNF, setParetoLogPDF, getParetoLogCDFNF, setParetoLogCDF, setParetoLogQuan, setParetoLogRand

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowetoLogPDFNF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPowetoLogPDFNF_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowetoLogPDFNF_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowetoLogPDFNF_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowetoLogPDFNF_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowetoLogPDFNF_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowetoLogPDFNF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowetoLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPowetoLogPDF_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowetoLogPDF_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowetoLogPDF_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowetoLogPDF_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowetoLogPDF_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowetoLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPowetoLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPowetoLogPDF_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPowetoLogPDF_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPowetoLogPDF_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPowetoLogPDF_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPowetoLogPDF_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPowetoLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowetoLogCDFNF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPowetoLogCDFNF_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowetoLogCDFNF_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowetoLogCDFNF_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowetoLogCDFNF_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowetoLogCDFNF_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowetoLogCDFNF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowetoLogCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPowetoLogCDF_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowetoLogCDF_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowetoLogCDF_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowetoLogCDF_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowetoLogCDF_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowetoLogCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPowetoLogCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPowetoLogCDF_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPowetoLogCDF_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPowetoLogCDF_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPowetoLogCDF_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPowetoLogCDF_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPowetoLogCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowetoLogQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPowetoLogQuan_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowetoLogQuan_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowetoLogQuan_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowetoLogQuan_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowetoLogQuan_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowetoLogQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPowetoLogQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPowetoLogQuan_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPowetoLogQuan_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPowetoLogQuan_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPowetoLogQuan_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPowetoLogQuan_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPowetoLogQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowetoLogRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPowetoLogRand_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowetoLogRand_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowetoLogRand_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowetoLogRand_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowetoLogRand_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowetoLogRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPowetoLogRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPowetoLogRand_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPowetoLogRand_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPowetoLogRand_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPowetoLogRand_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPowetoLogRand_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPowetoLogRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines
