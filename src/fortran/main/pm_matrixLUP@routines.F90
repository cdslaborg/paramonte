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
!>  This file contains procedure implementations of [pm_matrixLUP](@ref pm_matrixLUP).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_matrixLUP) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_lapack, only: lapackGETRF
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatLUP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IMP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SQM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatLUP_IMP_SQM_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatLUP_IMP_SQM_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatLUP_IMP_SQM_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatLUP_IMP_SQM_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixLUP@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatLUP_IMP_SQM_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixLUP@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatLUP_IMP_SQM_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatLUP_IMP_SQM_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatLUP_IMP_SQM_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatLUP_IMP_SQM_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixLUP@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatLUP_IMP_SQM_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixLUP@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SQM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
#define ITE_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatLUP_IMP_ITE_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatLUP_IMP_ITE_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatLUP_IMP_ITE_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatLUP_IMP_ITE_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatLUP_IMP_ITE_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ITE_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define REC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatLUP_IMP_REC_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatLUP_IMP_REC_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatLUP_IMP_REC_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatLUP_IMP_REC_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatLUP_IMP_REC_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef REC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IMP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define EXP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SQM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatLUP_EXP_SQM_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatLUP_EXP_SQM_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatLUP_EXP_SQM_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatLUP_EXP_SQM_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatLUP_EXP_SQM_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SQM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ITE_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatLUP_EXP_ITE_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatLUP_EXP_ITE_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatLUP_EXP_ITE_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatLUP_EXP_ITE_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatLUP_EXP_ITE_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ITE_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define REC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatLUP_EXP_REC_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatLUP_EXP_REC_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatLUP_EXP_REC_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatLUP_EXP_REC_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatLUP_EXP_REC_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixLUP@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef REC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef EXP_ENABLED
#endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatLUP_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines