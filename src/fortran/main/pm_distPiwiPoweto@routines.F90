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
!>  This file contains procedure implementations of [pm_distPiwiPoweto](@ref pm_distPiwiPoweto).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distPiwiPoweto) routines ! LCOV_EXCL_LINE

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

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPiwiPowetoLogPDFNF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPiwiPowetoLogPDFNFALD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPiwiPowetoLogPDFNFALD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPiwiPowetoLogPDFNFALD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPiwiPowetoLogPDFNFALD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPiwiPowetoLogPDFNFALD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPiwiPowetoLogPDFNFALC_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPiwiPowetoLogPDFNFALC_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPiwiPowetoLogPDFNFALC_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPiwiPowetoLogPDFNFALC_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPiwiPowetoLogPDFNFALC_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPiwiPowetoLogPDFNF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPiwiPowetoLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPiwiPowetoLogPDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPiwiPowetoLogPDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPiwiPowetoLogPDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPiwiPowetoLogPDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPiwiPowetoLogPDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPiwiPowetoLogPDF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPiwiPowetoLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPiwiPowetoLogPDFALL_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPiwiPowetoLogPDFALL_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPiwiPowetoLogPDFALL_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPiwiPowetoLogPDFALL_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPiwiPowetoLogPDFALL_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BAN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPiwiPowetoLogPDFBAN_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPiwiPowetoLogPDFBAN_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPiwiPowetoLogPDFBAN_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPiwiPowetoLogPDFBAN_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPiwiPowetoLogPDFBAN_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef BAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPiwiPowetoLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPiwiPowetoCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPiwiPowetoCDFALDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPiwiPowetoCDFALDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPiwiPowetoCDFALDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPiwiPowetoCDFALDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPiwiPowetoCDFALDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALLC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPiwiPowetoCDFALLC_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPiwiPowetoCDFALLC_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPiwiPowetoCDFALLC_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPiwiPowetoCDFALLC_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPiwiPowetoCDFALLC_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALLC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPiwiPowetoCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPiwiPowetoCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MAN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPiwiPowetoCDFMAN_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPiwiPowetoCDFMAN_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPiwiPowetoCDFMAN_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPiwiPowetoCDFMAN_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPiwiPowetoCDFMAN_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BAN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPiwiPowetoCDFBAN_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPiwiPowetoCDFBAN_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPiwiPowetoCDFBAN_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPiwiPowetoCDFBAN_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPiwiPowetoCDFBAN_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef BAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPiwiPowetoCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines
