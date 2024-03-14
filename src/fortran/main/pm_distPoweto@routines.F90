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

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowetoLogPDFNF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPowetoLogPDFNFALD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowetoLogPDFNFALD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowetoLogPDFNFALD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowetoLogPDFNFALD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowetoLogPDFNFALD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
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
    module procedure getPowetoLogPDFNFALC_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowetoLogPDFNFALC_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowetoLogPDFNFALC_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowetoLogPDFNFALC_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowetoLogPDFNFALC_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowetoLogPDFNF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowetoLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPowetoLogPDF_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowetoLogPDF_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowetoLogPDF_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowetoLogPDF_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowetoLogPDF_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowetoLogPDF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPowetoLogPDF_ENABLED 1

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
    module procedure setPowetoLogPDFALL_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPowetoLogPDFALL_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPowetoLogPDFALL_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPowetoLogPDFALL_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPowetoLogPDFALL_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
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
    module procedure setPowetoLogPDFBAN_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPowetoLogPDFBAN_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPowetoLogPDFBAN_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPowetoLogPDFBAN_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPowetoLogPDFBAN_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
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

#undef setPowetoLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowetoCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPowetoCDFALDD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowetoCDFALDD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowetoCDFALDD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowetoCDFALDD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowetoCDFALDD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
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
    module procedure getPowetoCDFALLC_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowetoCDFALLC_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowetoCDFALLC_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowetoCDFALLC_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowetoCDFALLC_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALLC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowetoCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPowetoCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MAN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPowetoCDFMAN_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPowetoCDFMAN_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPowetoCDFMAN_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPowetoCDFMAN_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPowetoCDFMAN_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
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
    module procedure setPowetoCDFBAN_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPowetoCDFBAN_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPowetoCDFBAN_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPowetoCDFBAN_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPowetoCDFBAN_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef BAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPowetoCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines
