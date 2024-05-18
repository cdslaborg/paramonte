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
!>  This file contains procedure implementations of [pm_distKolm](@ref pm_distKolm).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distKolm) routines ! LCOV_EXCL_LINE

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

#define getKolmPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getKolmPDF_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getKolmPDF_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getKolmPDF_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getKolmPDF_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getKolmPDF_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getKolmPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setKolmPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setKolmPDF_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setKolmPDF_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setKolmPDF_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setKolmPDF_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setKolmPDF_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setKolmPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getKolmCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getKolmCDF_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getKolmCDF_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getKolmCDF_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getKolmCDF_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getKolmCDF_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getKolmCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setKolmCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setKolmCDF_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setKolmCDF_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setKolmCDF_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setKolmCDF_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setKolmCDF_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setKolmCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getKolmQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getKolmQuan_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getKolmQuan_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getKolmQuan_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getKolmQuan_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getKolmQuan_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getKolmQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setKolmQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setKolmQuan_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setKolmQuan_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setKolmQuan_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setKolmQuan_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setKolmQuan_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setKolmQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getKolmRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getKolmRand_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getKolmRand_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getKolmRand_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getKolmRand_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getKolmRand_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getKolmRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setKolmRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setKolmRand_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setKolmRand_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setKolmRand_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setKolmRand_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setKolmRand_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setKolmRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines