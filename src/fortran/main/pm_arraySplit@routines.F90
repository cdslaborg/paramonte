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
!>  This file contains procedure implementations of [pm_arraySplit](@ref pm_arraySplit).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_arraySplit) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_arrayUnique, only: getUnique
    use pm_arraySort, only: setSorted

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define setSplit_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define Fixed_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define DefCom_ENABLED 1
!#define DefIns_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define D0_D0_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setSplitFixDefComDefIns_D0_D0_SK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setSplitFixDefComDefIns_D0_D0_SK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setSplitFixDefComDefIns_D0_D0_SK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setSplitFixDefComDefIns_D0_D0_SK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setSplitFixDefComDefIns_D0_D0_SK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef D0_D0_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define D1_D0_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_SK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_SK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_SK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_SK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_SK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_IK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_IK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_IK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_IK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_IK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define LK_ENABLED 1
!
!#if LK5_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_LK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK4_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_LK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK3_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_LK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK2_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_LK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK1_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_LK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef LK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define CK_ENABLED 1
!
!#if CK5_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_CK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK4_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_CK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK3_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_CK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK2_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_CK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK1_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_CK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef CK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_RK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_RK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_RK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_RK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D0_RK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef D1_D0_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define D1_D1_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_SK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_SK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_SK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_SK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_SK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_IK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_IK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_IK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_IK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_IK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define LK_ENABLED 1
!
!#if LK5_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_LK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK4_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_LK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK3_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_LK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK2_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_LK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK1_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_LK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef LK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define CK_ENABLED 1
!
!#if CK5_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_CK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK4_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_CK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK3_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_CK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK2_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_CK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK1_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_CK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef CK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_RK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_RK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_RK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_RK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure setSplitFixDefComDefIns_D1_D1_RK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef D1_D1_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef DefCom_ENABLED
!#undef DefIns_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if FORTRAN_ENABLED
!#define CusCom_ENABLED 1
!#define DefIns_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define D0_D0_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setSplitFixCusComDefIns_D0_D0_SK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setSplitFixCusComDefIns_D0_D0_SK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setSplitFixCusComDefIns_D0_D0_SK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setSplitFixCusComDefIns_D0_D0_SK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setSplitFixCusComDefIns_D0_D0_SK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef D0_D0_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define D1_D0_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_SK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_SK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_SK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_SK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_SK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_IK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_IK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_IK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_IK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_IK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define LK_ENABLED 1
!
!#if LK5_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_LK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK4_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_LK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK3_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_LK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK2_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_LK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK1_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_LK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef LK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define CK_ENABLED 1
!
!#if CK5_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_CK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK4_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_CK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK3_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_CK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK2_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_CK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK1_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_CK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef CK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_RK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_RK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_RK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_RK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D0_RK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef D1_D0_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define D1_D1_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_SK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_SK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_SK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_SK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_SK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_IK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_IK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_IK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_IK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_IK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define LK_ENABLED 1
!
!#if LK5_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_LK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK4_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_LK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK3_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_LK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK2_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_LK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK1_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_LK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef LK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define CK_ENABLED 1
!
!#if CK5_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_CK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK4_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_CK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK3_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_CK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK2_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_CK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK1_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_CK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef CK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_RK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_RK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_RK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_RK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure setSplitFixCusComDefIns_D1_D1_RK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef D1_D1_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef CusCom_ENABLED
!#endif
!FORTRAN_ENABLED
!#undef DefIns_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define DefCom_ENABLED 1
!#define CusIns_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define D0_D0_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setSplitFixDefComCusIns_D0_D0_SK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setSplitFixDefComCusIns_D0_D0_SK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setSplitFixDefComCusIns_D0_D0_SK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setSplitFixDefComCusIns_D0_D0_SK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setSplitFixDefComCusIns_D0_D0_SK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef D0_D0_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define D1_D0_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_SK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_SK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_SK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_SK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_SK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_IK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_IK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_IK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_IK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_IK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define LK_ENABLED 1
!
!#if LK5_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_LK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK4_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_LK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK3_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_LK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK2_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_LK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK1_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_LK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef LK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define CK_ENABLED 1
!
!#if CK5_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_CK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK4_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_CK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK3_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_CK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK2_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_CK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK1_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_CK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef CK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_RK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_RK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_RK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_RK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D0_RK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef D1_D0_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define D1_D1_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_SK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_SK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_SK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_SK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_SK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_IK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_IK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_IK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_IK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_IK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define LK_ENABLED 1
!
!#if LK5_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_LK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK4_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_LK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK3_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_LK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK2_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_LK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK1_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_LK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef LK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define CK_ENABLED 1
!
!#if CK5_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_CK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK4_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_CK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK3_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_CK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK2_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_CK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK1_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_CK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef CK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_RK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_RK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_RK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_RK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure setSplitFixDefComCusIns_D1_D1_RK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef D1_D1_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef CusIns_ENABLED
!#undef DefCom_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if FORTRAN_ENABLED
!#define CusCom_ENABLED 1
!#define CusIns_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define D0_D0_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setSplitFixCusComCusIns_D0_D0_SK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setSplitFixCusComCusIns_D0_D0_SK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setSplitFixCusComCusIns_D0_D0_SK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setSplitFixCusComCusIns_D0_D0_SK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setSplitFixCusComCusIns_D0_D0_SK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef D0_D0_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define D1_D0_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_SK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_SK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_SK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_SK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_SK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_IK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_IK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_IK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_IK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_IK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define LK_ENABLED 1
!
!#if LK5_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_LK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK4_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_LK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK3_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_LK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK2_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_LK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK1_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_LK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef LK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define CK_ENABLED 1
!
!#if CK5_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_CK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK4_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_CK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK3_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_CK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK2_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_CK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK1_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_CK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef CK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_RK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_RK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_RK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_RK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D0_RK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef D1_D0_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define D1_D1_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_SK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_SK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_SK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_SK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_SK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_IK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_IK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_IK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_IK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_IK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define LK_ENABLED 1
!
!#if LK5_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_LK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK4_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_LK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK3_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_LK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK2_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_LK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if LK1_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_LK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef LK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define CK_ENABLED 1
!
!#if CK5_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_CK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK4_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_CK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK3_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_CK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK2_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_CK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK1_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_CK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef CK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_RK5
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_RK4
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_RK3
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_RK2
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure setSplitFixCusComCusIns_D1_D1_RK1
!#include "pm_arraySplit@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef D1_D1_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef CusIns_ENABLED
!#undef CusCom_ENABLED
!#endif
!FORTRAN_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef Fixed_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef setSplit_ENABLED
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSplit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Index_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1
#define DefIns_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitIndDefComDefIns_D0_D0_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitIndDefComDefIns_D0_D0_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitIndDefComDefIns_D0_D0_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitIndDefComDefIns_D0_D0_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitIndDefComDefIns_D0_D0_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_IK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_IK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_IK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_IK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_IK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_LK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_LK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_LK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_LK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_LK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_CK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_CK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_CK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_CK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_CK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_RK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_RK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_RK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_RK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D0_RK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_IK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_IK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_IK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_IK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_IK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_LK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_LK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_LK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_LK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_LK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_CK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_CK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_CK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_CK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_CK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_RK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_RK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_RK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_RK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSplitIndDefComDefIns_D1_D1_RK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED
#undef DefIns_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1
#define DefIns_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitIndCusComDefIns_D0_D0_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitIndCusComDefIns_D0_D0_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitIndCusComDefIns_D0_D0_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitIndCusComDefIns_D0_D0_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitIndCusComDefIns_D0_D0_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_IK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_IK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_IK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_IK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_IK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_LK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_LK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_LK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_LK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_LK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_CK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_CK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_CK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_CK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_CK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_RK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_RK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_RK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_RK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D0_RK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_IK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_IK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_IK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_IK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_IK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_LK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_LK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_LK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_LK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_LK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_CK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_CK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_CK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_CK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_CK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_RK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_RK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_RK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_RK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSplitIndCusComDefIns_D1_D1_RK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED
#undef DefIns_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1
#if FORTRAN_ENABLED
#define CusIns_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitIndDefComCusIns_D0_D0_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitIndDefComCusIns_D0_D0_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitIndDefComCusIns_D0_D0_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitIndDefComCusIns_D0_D0_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitIndDefComCusIns_D0_D0_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_IK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_IK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_IK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_IK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_IK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_LK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_LK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_LK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_LK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_LK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_CK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_CK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_CK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_CK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_CK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_RK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_RK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_RK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_RK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D0_RK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_IK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_IK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_IK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_IK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_IK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_LK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_LK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_LK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_LK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_LK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_CK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_CK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_CK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_CK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_CK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_RK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_RK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_RK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_RK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSplitIndDefComCusIns_D1_D1_RK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusIns_ENABLED
#endif
!FORTRAN_ENABLED
#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1
#if FORTRAN_ENABLED
#define CusIns_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitIndCusComCusIns_D0_D0_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitIndCusComCusIns_D0_D0_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitIndCusComCusIns_D0_D0_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitIndCusComCusIns_D0_D0_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitIndCusComCusIns_D0_D0_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_IK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_IK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_IK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_IK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_IK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_LK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_LK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_LK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_LK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_LK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_CK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_CK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_CK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_CK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_CK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_RK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_RK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_RK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_RK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D0_RK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_IK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_IK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_IK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_IK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_IK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_LK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_LK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_LK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_LK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_LK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_CK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_CK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_CK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_CK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_CK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_RK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_RK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_RK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_RK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSplitIndCusComCusIns_D1_D1_RK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusIns_ENABLED
#endif
!FORTRAN_ENABLED
#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Index_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSplit_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSplit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define Jagged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefIns_ENABLED 1
#define DefCom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitConDefComDefIns_D0_D0_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitConDefComDefIns_D0_D0_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitConDefComDefIns_D0_D0_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitConDefComDefIns_D0_D0_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitConDefComDefIns_D0_D0_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_IK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_IK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_IK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_IK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_IK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_LK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_LK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_LK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_LK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_LK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_CK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_CK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_CK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_CK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_CK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_RK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_RK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_RK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_RK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSplitConDefComDefIns_D1_D0_RK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_IK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_IK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_IK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_IK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_IK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_LK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_LK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_LK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_LK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_LK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_CK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_CK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_CK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_CK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_CK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_RK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_RK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_RK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_RK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSplitConDefComDefIns_D1_D1_RK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefIns_ENABLED
#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1
#define DefIns_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitConCusComDefIns_D0_D0_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitConCusComDefIns_D0_D0_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitConCusComDefIns_D0_D0_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitConCusComDefIns_D0_D0_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitConCusComDefIns_D0_D0_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_IK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_IK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_IK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_IK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_IK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_LK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_LK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_LK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_LK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_LK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_CK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_CK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_CK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_CK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_CK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_RK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_RK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_RK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_RK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSplitConCusComDefIns_D1_D0_RK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_IK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_IK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_IK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_IK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_IK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_LK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_LK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_LK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_LK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_LK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_CK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_CK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_CK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_CK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_CK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_RK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_RK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_RK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_RK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSplitConCusComDefIns_D1_D1_RK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED
#undef DefIns_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1
#if FORTRAN_ENABLED
#define CusIns_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitConDefComCusIns_D0_D0_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitConDefComCusIns_D0_D0_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitConDefComCusIns_D0_D0_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitConDefComCusIns_D0_D0_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitConDefComCusIns_D0_D0_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_IK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_IK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_IK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_IK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_IK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_LK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_LK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_LK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_LK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_LK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_CK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_CK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_CK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_CK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_CK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_RK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_RK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_RK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_RK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSplitConDefComCusIns_D1_D0_RK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_IK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_IK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_IK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_IK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_IK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_LK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_LK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_LK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_LK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_LK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_CK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_CK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_CK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_CK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_CK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_RK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_RK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_RK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_RK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSplitConDefComCusIns_D1_D1_RK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusIns_ENABLED
#endif
!FORTRAN_ENABLED
#undef DefCom_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1
#if FORTRAN_ENABLED
#define CusIns_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitConCusComCusIns_D0_D0_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitConCusComCusIns_D0_D0_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitConCusComCusIns_D0_D0_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitConCusComCusIns_D0_D0_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitConCusComCusIns_D0_D0_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_IK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_IK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_IK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_IK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_IK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_LK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_LK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_LK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_LK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_LK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_CK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_CK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_CK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_CK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_CK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_RK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_RK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_RK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_RK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSplitConCusComCusIns_D1_D0_RK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_SK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_SK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_SK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_SK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_SK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_IK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_IK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_IK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_IK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_IK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_LK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_LK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_LK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_LK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_LK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_CK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_CK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_CK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_CK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_CK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_RK5
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_RK4
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_RK3
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_RK2
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSplitConCusComCusIns_D1_D1_RK1
#include "pm_arraySplit@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusIns_ENABLED
#endif
!FORTRAN_ENABLED
#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Jagged_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSplit_ENABLED 1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSplit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Jagged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1
#define DefIns_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure setSplitBoxDefComDefIns_D0_D0_SK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure setSplitBoxDefComDefIns_D1_D0_SK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure setSplitBoxDefComDefIns_D1_D0_IK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure setSplitBoxDefComDefIns_D1_D0_LK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure setSplitBoxDefComDefIns_D1_D0_CK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure setSplitBoxDefComDefIns_D1_D0_RK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure setSplitBoxDefComDefIns_D1_D1_SK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure setSplitBoxDefComDefIns_D1_D1_IK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure setSplitBoxDefComDefIns_D1_D1_LK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure setSplitBoxDefComDefIns_D1_D1_CK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure setSplitBoxDefComDefIns_D1_D1_RK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED
#undef DefIns_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1
#define DefIns_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure setSplitBoxCusComDefIns_D0_D0_SK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure setSplitBoxCusComDefIns_D1_D0_SK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure setSplitBoxCusComDefIns_D1_D0_IK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure setSplitBoxCusComDefIns_D1_D0_LK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure setSplitBoxCusComDefIns_D1_D0_CK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure setSplitBoxCusComDefIns_D1_D0_RK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure setSplitBoxCusComDefIns_D1_D1_SK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure setSplitBoxCusComDefIns_D1_D1_IK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure setSplitBoxCusComDefIns_D1_D1_LK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure setSplitBoxCusComDefIns_D1_D1_CK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure setSplitBoxCusComDefIns_D1_D1_RK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED
#undef DefIns_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DefCom_ENABLED 1
#if FORTRAN_ENABLED
#define CusIns_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure setSplitBoxDefComCusIns_D0_D0_SK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure setSplitBoxDefComCusIns_D1_D0_SK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure setSplitBoxDefComCusIns_D1_D0_IK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure setSplitBoxDefComCusIns_D1_D0_LK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure setSplitBoxDefComCusIns_D1_D0_CK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure setSplitBoxDefComCusIns_D1_D0_RK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure setSplitBoxDefComCusIns_D1_D1_SK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure setSplitBoxDefComCusIns_D1_D1_IK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure setSplitBoxDefComCusIns_D1_D1_LK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure setSplitBoxDefComCusIns_D1_D1_CK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure setSplitBoxDefComCusIns_D1_D1_RK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DefCom_ENABLED
#undef CusIns_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define CusCom_ENABLED 1
#if FORTRAN_ENABLED
#define CusIns_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure setSplitBoxCusComCusIns_D0_D0_SK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure setSplitBoxCusComCusIns_D1_D0_SK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure setSplitBoxCusComCusIns_D1_D0_IK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure setSplitBoxCusComCusIns_D1_D0_LK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure setSplitBoxCusComCusIns_D1_D0_CK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure setSplitBoxCusComCusIns_D1_D0_RK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure setSplitBoxCusComCusIns_D1_D1_SK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure setSplitBoxCusComCusIns_D1_D1_IK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure setSplitBoxCusComCusIns_D1_D1_LK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure setSplitBoxCusComCusIns_D1_D1_CK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure setSplitBoxCusComCusIns_D1_D1_RK
#include "pm_arraySplit@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CusIns_ENABLED
#endif
!FORTRAN_ENABLED
#undef CusCom_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Jagged_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSplit_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines