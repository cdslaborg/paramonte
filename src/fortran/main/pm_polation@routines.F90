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
!>  This file contains procedure implementations of [pm_polation](@ref pm_polation).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_polation) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_arraySort, only: isAscending, isAscendingAll
    use pm_arraySearch, only: getBin
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getExtrap_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ND1_ENABLED 1

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
    module procedure getExtrapPWLN_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExtrapPWLN_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExtrapPWLN_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExtrapPWLN_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExtrapPWLN_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getExtrapPWLN_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExtrapPWLN_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExtrapPWLN_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExtrapPWLN_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExtrapPWLN_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getExtrapMEAN_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExtrapMEAN_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExtrapMEAN_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExtrapMEAN_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExtrapMEAN_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getExtrapMEAN_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExtrapMEAN_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExtrapMEAN_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExtrapMEAN_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExtrapMEAN_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getExtrapNEAR_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExtrapNEAR_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExtrapNEAR_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExtrapNEAR_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExtrapNEAR_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getExtrapNEAR_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExtrapNEAR_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExtrapNEAR_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExtrapNEAR_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExtrapNEAR_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getExtrapNEXT_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExtrapNEXT_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExtrapNEXT_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExtrapNEXT_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExtrapNEXT_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getExtrapNEXT_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExtrapNEXT_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExtrapNEXT_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExtrapNEXT_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExtrapNEXT_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getExtrapPREV_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExtrapPREV_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExtrapPREV_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExtrapPREV_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExtrapPREV_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getExtrapPREV_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExtrapPREV_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExtrapPREV_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExtrapPREV_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExtrapPREV_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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

#define MNPLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getExtrapMNPLD_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExtrapMNPLD_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExtrapMNPLD_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExtrapMNPLD_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExtrapMNPLD_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getExtrapMNPLD_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExtrapMNPLD_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExtrapMNPLD_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExtrapMNPLD_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExtrapMNPLD_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MNPLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ND1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getExtrap_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setExtrap_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ND1_ENABLED 1

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
    module procedure setExtrapPWLN_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExtrapPWLN_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExtrapPWLN_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExtrapPWLN_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExtrapPWLN_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setExtrapPWLN_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExtrapPWLN_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExtrapPWLN_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExtrapPWLN_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExtrapPWLN_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setExtrapMEAN_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExtrapMEAN_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExtrapMEAN_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExtrapMEAN_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExtrapMEAN_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setExtrapMEAN_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExtrapMEAN_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExtrapMEAN_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExtrapMEAN_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExtrapMEAN_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setExtrapNEAR_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExtrapNEAR_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExtrapNEAR_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExtrapNEAR_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExtrapNEAR_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setExtrapNEAR_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExtrapNEAR_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExtrapNEAR_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExtrapNEAR_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExtrapNEAR_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setExtrapNEXT_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExtrapNEXT_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExtrapNEXT_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExtrapNEXT_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExtrapNEXT_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setExtrapNEXT_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExtrapNEXT_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExtrapNEXT_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExtrapNEXT_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExtrapNEXT_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setExtrapPREV_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExtrapPREV_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExtrapPREV_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExtrapPREV_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExtrapPREV_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setExtrapPREV_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExtrapPREV_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExtrapPREV_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExtrapPREV_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExtrapPREV_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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

#define MNPLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setExtrapMNPLD_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExtrapMNPLD_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExtrapMNPLD_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExtrapMNPLD_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExtrapMNPLD_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setExtrapMNPLD_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExtrapMNPLD_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExtrapMNPLD_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExtrapMNPLD_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExtrapMNPLD_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MNPLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MNPLE_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setExtrapMNPLE_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExtrapMNPLE_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExtrapMNPLE_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExtrapMNPLE_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExtrapMNPLE_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setExtrapMNPLE_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExtrapMNPLE_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExtrapMNPLE_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExtrapMNPLE_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExtrapMNPLE_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MNPLE_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ND1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setExtrap_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getInterp_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ND1_ENABLED 1

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
    module procedure getInterpPWLN_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getInterpPWLN_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getInterpPWLN_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getInterpPWLN_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getInterpPWLN_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getInterpPWLN_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getInterpPWLN_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getInterpPWLN_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getInterpPWLN_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getInterpPWLN_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getInterpMEAN_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getInterpMEAN_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getInterpMEAN_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getInterpMEAN_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getInterpMEAN_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getInterpMEAN_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getInterpMEAN_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getInterpMEAN_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getInterpMEAN_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getInterpMEAN_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getInterpNEAR_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getInterpNEAR_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getInterpNEAR_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getInterpNEAR_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getInterpNEAR_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getInterpNEAR_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getInterpNEAR_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getInterpNEAR_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getInterpNEAR_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getInterpNEAR_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getInterpNEXT_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getInterpNEXT_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getInterpNEXT_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getInterpNEXT_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getInterpNEXT_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getInterpNEXT_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getInterpNEXT_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getInterpNEXT_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getInterpNEXT_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getInterpNEXT_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getInterpPREV_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getInterpPREV_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getInterpPREV_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getInterpPREV_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getInterpPREV_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getInterpPREV_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getInterpPREV_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getInterpPREV_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getInterpPREV_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getInterpPREV_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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

#define MNPLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getInterpMNPLD_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getInterpMNPLD_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getInterpMNPLD_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getInterpMNPLD_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getInterpMNPLD_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure getInterpMNPLD_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getInterpMNPLD_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getInterpMNPLD_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getInterpMNPLD_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getInterpMNPLD_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MNPLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ND1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getInterp_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setInterp_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ND1_ENABLED 1

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
    module procedure setInterpPWLN_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setInterpPWLN_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setInterpPWLN_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setInterpPWLN_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setInterpPWLN_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setInterpPWLN_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setInterpPWLN_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setInterpPWLN_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setInterpPWLN_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setInterpPWLN_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setInterpMEAN_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setInterpMEAN_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setInterpMEAN_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setInterpMEAN_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setInterpMEAN_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setInterpMEAN_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setInterpMEAN_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setInterpMEAN_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setInterpMEAN_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setInterpMEAN_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setInterpNEAR_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setInterpNEAR_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setInterpNEAR_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setInterpNEAR_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setInterpNEAR_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setInterpNEAR_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setInterpNEAR_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setInterpNEAR_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setInterpNEAR_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setInterpNEAR_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setInterpNEXT_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setInterpNEXT_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setInterpNEXT_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setInterpNEXT_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setInterpNEXT_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setInterpNEXT_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setInterpNEXT_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setInterpNEXT_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setInterpNEXT_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setInterpNEXT_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setInterpPREV_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setInterpPREV_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setInterpPREV_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setInterpPREV_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setInterpPREV_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setInterpPREV_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setInterpPREV_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setInterpPREV_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setInterpPREV_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setInterpPREV_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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

#define MNPLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setInterpMNPLD_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setInterpMNPLD_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setInterpMNPLD_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setInterpMNPLD_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setInterpMNPLD_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setInterpMNPLD_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setInterpMNPLD_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setInterpMNPLD_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setInterpMNPLD_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setInterpMNPLD_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MNPLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MNPLE_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QD0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setInterpMNPLE_ND1_QD0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setInterpMNPLE_ND1_QD0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setInterpMNPLE_ND1_QD0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setInterpMNPLE_ND1_QD0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setInterpMNPLE_ND1_QD0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
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
    module procedure setInterpMNPLE_ND1_QD1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setInterpMNPLE_ND1_QD1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setInterpMNPLE_ND1_QD1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setInterpMNPLE_ND1_QD1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setInterpMNPLE_ND1_QD1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_polation@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MNPLE_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ND1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setInterp_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines