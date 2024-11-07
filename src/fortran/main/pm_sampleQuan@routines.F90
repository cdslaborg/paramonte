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
!>  This file contains procedure implementations of [pm_sampleQuan](@ref pm_sampleQuan).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_sampleQuan) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_arraySort, only: setSorted
    use pm_sampleECDF, only: setECDF
    use pm_polation, only: setExtrap
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ND1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PWLN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPWLN_ND1_QD0_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPWLN_ND1_QD0_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPWLN_ND1_QD0_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPWLN_ND1_QD0_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPWLN_ND1_QD0_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPWLN_ND1_QD1_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPWLN_ND1_QD1_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPWLN_ND1_QD1_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPWLN_ND1_QD1_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPWLN_ND1_QD1_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PWLN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MEAN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanMEAN_ND1_QD0_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanMEAN_ND1_QD0_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanMEAN_ND1_QD0_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanMEAN_ND1_QD0_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanMEAN_ND1_QD0_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanMEAN_ND1_QD1_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanMEAN_ND1_QD1_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanMEAN_ND1_QD1_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanMEAN_ND1_QD1_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanMEAN_ND1_QD1_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MEAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NEAR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEAR_ND1_QD0_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEAR_ND1_QD0_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEAR_ND1_QD0_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEAR_ND1_QD0_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEAR_ND1_QD0_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEAR_ND1_QD1_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEAR_ND1_QD1_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEAR_ND1_QD1_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEAR_ND1_QD1_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEAR_ND1_QD1_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NEAR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NEXT_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEXT_ND1_QD0_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEXT_ND1_QD0_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEXT_ND1_QD0_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEXT_ND1_QD0_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEXT_ND1_QD0_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEXT_ND1_QD1_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEXT_ND1_QD1_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEXT_ND1_QD1_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEXT_ND1_QD1_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEXT_ND1_QD1_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NEXT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PREV_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPREV_ND1_QD0_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPREV_ND1_QD0_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPREV_ND1_QD0_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPREV_ND1_QD0_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPREV_ND1_QD0_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPREV_ND1_QD1_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPREV_ND1_QD1_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPREV_ND1_QD1_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPREV_ND1_QD1_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPREV_ND1_QD1_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PREV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PWLN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPWLN_ND1_QD0_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPWLN_ND1_QD0_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPWLN_ND1_QD0_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPWLN_ND1_QD0_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPWLN_ND1_QD0_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPWLN_ND1_QD1_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPWLN_ND1_QD1_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPWLN_ND1_QD1_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPWLN_ND1_QD1_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPWLN_ND1_QD1_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PWLN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MEAN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanMEAN_ND1_QD0_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanMEAN_ND1_QD0_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanMEAN_ND1_QD0_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanMEAN_ND1_QD0_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanMEAN_ND1_QD0_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanMEAN_ND1_QD1_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanMEAN_ND1_QD1_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanMEAN_ND1_QD1_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanMEAN_ND1_QD1_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanMEAN_ND1_QD1_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MEAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NEAR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEAR_ND1_QD0_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEAR_ND1_QD0_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEAR_ND1_QD0_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEAR_ND1_QD0_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEAR_ND1_QD0_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEAR_ND1_QD1_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEAR_ND1_QD1_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEAR_ND1_QD1_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEAR_ND1_QD1_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEAR_ND1_QD1_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NEAR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NEXT_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEXT_ND1_QD0_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEXT_ND1_QD0_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEXT_ND1_QD0_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEXT_ND1_QD0_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEXT_ND1_QD0_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEXT_ND1_QD1_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEXT_ND1_QD1_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEXT_ND1_QD1_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEXT_ND1_QD1_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEXT_ND1_QD1_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NEXT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PREV_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPREV_ND1_QD0_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPREV_ND1_QD0_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPREV_ND1_QD0_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPREV_ND1_QD0_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPREV_ND1_QD0_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPREV_ND1_QD1_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPREV_ND1_QD1_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPREV_ND1_QD1_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPREV_ND1_QD1_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPREV_ND1_QD1_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PREV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PWLN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPWLN_ND1_QD0_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPWLN_ND1_QD0_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPWLN_ND1_QD0_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPWLN_ND1_QD0_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPWLN_ND1_QD0_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPWLN_ND1_QD1_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPWLN_ND1_QD1_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPWLN_ND1_QD1_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPWLN_ND1_QD1_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPWLN_ND1_QD1_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PWLN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MEAN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanMEAN_ND1_QD0_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanMEAN_ND1_QD0_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanMEAN_ND1_QD0_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanMEAN_ND1_QD0_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanMEAN_ND1_QD0_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanMEAN_ND1_QD1_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanMEAN_ND1_QD1_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanMEAN_ND1_QD1_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanMEAN_ND1_QD1_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanMEAN_ND1_QD1_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MEAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NEAR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEAR_ND1_QD0_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEAR_ND1_QD0_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEAR_ND1_QD0_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEAR_ND1_QD0_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEAR_ND1_QD0_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEAR_ND1_QD1_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEAR_ND1_QD1_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEAR_ND1_QD1_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEAR_ND1_QD1_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEAR_ND1_QD1_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NEAR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NEXT_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEXT_ND1_QD0_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEXT_ND1_QD0_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEXT_ND1_QD0_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEXT_ND1_QD0_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEXT_ND1_QD0_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEXT_ND1_QD1_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEXT_ND1_QD1_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEXT_ND1_QD1_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEXT_ND1_QD1_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEXT_ND1_QD1_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NEXT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PREV_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPREV_ND1_QD0_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPREV_ND1_QD0_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPREV_ND1_QD0_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPREV_ND1_QD0_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPREV_ND1_QD0_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPREV_ND1_QD1_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPREV_ND1_QD1_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPREV_ND1_QD1_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPREV_ND1_QD1_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPREV_ND1_QD1_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PREV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ND1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ND2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PWLN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPWLN_ND2_QD0_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPWLN_ND2_QD0_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPWLN_ND2_QD0_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPWLN_ND2_QD0_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPWLN_ND2_QD0_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPWLN_ND2_QD1_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPWLN_ND2_QD1_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPWLN_ND2_QD1_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPWLN_ND2_QD1_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPWLN_ND2_QD1_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PWLN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MEAN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanMEAN_ND2_QD0_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanMEAN_ND2_QD0_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanMEAN_ND2_QD0_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanMEAN_ND2_QD0_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanMEAN_ND2_QD0_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanMEAN_ND2_QD1_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanMEAN_ND2_QD1_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanMEAN_ND2_QD1_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanMEAN_ND2_QD1_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanMEAN_ND2_QD1_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MEAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NEAR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEAR_ND2_QD0_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEAR_ND2_QD0_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEAR_ND2_QD0_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEAR_ND2_QD0_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEAR_ND2_QD0_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEAR_ND2_QD1_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEAR_ND2_QD1_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEAR_ND2_QD1_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEAR_ND2_QD1_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEAR_ND2_QD1_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NEAR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NEXT_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEXT_ND2_QD0_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEXT_ND2_QD0_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEXT_ND2_QD0_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEXT_ND2_QD0_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEXT_ND2_QD0_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEXT_ND2_QD1_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEXT_ND2_QD1_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEXT_ND2_QD1_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEXT_ND2_QD1_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEXT_ND2_QD1_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NEXT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PREV_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPREV_ND2_QD0_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPREV_ND2_QD0_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPREV_ND2_QD0_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPREV_ND2_QD0_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPREV_ND2_QD0_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPREV_ND2_QD1_WNO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPREV_ND2_QD1_WNO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPREV_ND2_QD1_WNO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPREV_ND2_QD1_WNO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPREV_ND2_QD1_WNO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PREV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PWLN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPWLN_ND2_QD0_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPWLN_ND2_QD0_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPWLN_ND2_QD0_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPWLN_ND2_QD0_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPWLN_ND2_QD0_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPWLN_ND2_QD1_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPWLN_ND2_QD1_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPWLN_ND2_QD1_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPWLN_ND2_QD1_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPWLN_ND2_QD1_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PWLN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MEAN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanMEAN_ND2_QD0_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanMEAN_ND2_QD0_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanMEAN_ND2_QD0_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanMEAN_ND2_QD0_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanMEAN_ND2_QD0_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanMEAN_ND2_QD1_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanMEAN_ND2_QD1_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanMEAN_ND2_QD1_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanMEAN_ND2_QD1_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanMEAN_ND2_QD1_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MEAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NEAR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEAR_ND2_QD0_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEAR_ND2_QD0_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEAR_ND2_QD0_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEAR_ND2_QD0_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEAR_ND2_QD0_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEAR_ND2_QD1_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEAR_ND2_QD1_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEAR_ND2_QD1_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEAR_ND2_QD1_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEAR_ND2_QD1_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NEAR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NEXT_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEXT_ND2_QD0_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEXT_ND2_QD0_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEXT_ND2_QD0_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEXT_ND2_QD0_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEXT_ND2_QD0_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEXT_ND2_QD1_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEXT_ND2_QD1_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEXT_ND2_QD1_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEXT_ND2_QD1_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEXT_ND2_QD1_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NEXT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PREV_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPREV_ND2_QD0_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPREV_ND2_QD0_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPREV_ND2_QD0_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPREV_ND2_QD0_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPREV_ND2_QD0_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPREV_ND2_QD1_WTI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPREV_ND2_QD1_WTI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPREV_ND2_QD1_WTI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPREV_ND2_QD1_WTI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPREV_ND2_QD1_WTI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PREV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PWLN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPWLN_ND2_QD0_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPWLN_ND2_QD0_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPWLN_ND2_QD0_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPWLN_ND2_QD0_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPWLN_ND2_QD0_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPWLN_ND2_QD1_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPWLN_ND2_QD1_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPWLN_ND2_QD1_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPWLN_ND2_QD1_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPWLN_ND2_QD1_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PWLN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MEAN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanMEAN_ND2_QD0_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanMEAN_ND2_QD0_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanMEAN_ND2_QD0_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanMEAN_ND2_QD0_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanMEAN_ND2_QD0_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanMEAN_ND2_QD1_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanMEAN_ND2_QD1_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanMEAN_ND2_QD1_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanMEAN_ND2_QD1_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanMEAN_ND2_QD1_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MEAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NEAR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEAR_ND2_QD0_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEAR_ND2_QD0_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEAR_ND2_QD0_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEAR_ND2_QD0_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEAR_ND2_QD0_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEAR_ND2_QD1_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEAR_ND2_QD1_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEAR_ND2_QD1_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEAR_ND2_QD1_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEAR_ND2_QD1_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NEAR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NEXT_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEXT_ND2_QD0_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEXT_ND2_QD0_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEXT_ND2_QD0_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEXT_ND2_QD0_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEXT_ND2_QD0_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanNEXT_ND2_QD1_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanNEXT_ND2_QD1_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanNEXT_ND2_QD1_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanNEXT_ND2_QD1_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanNEXT_ND2_QD1_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NEXT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PREV_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPREV_ND2_QD0_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPREV_ND2_QD0_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPREV_ND2_QD0_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPREV_ND2_QD0_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPREV_ND2_QD0_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuanPREV_ND2_QD1_WTR_RK5
        use pm_kind, only: TKG => RK5
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuanPREV_ND2_QD1_WTR_RK4
        use pm_kind, only: TKG => RK4
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuanPREV_ND2_QD1_WTR_RK3
        use pm_kind, only: TKG => RK3
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuanPREV_ND2_QD1_WTR_RK2
        use pm_kind, only: TKG => RK2
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuanPREV_ND2_QD1_WTR_RK1
        use pm_kind, only: TKG => RK1
#include "pm_sampleQuan@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PREV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ND2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines