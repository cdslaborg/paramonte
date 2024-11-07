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
!>  This file contains procedure implementations of [pm_matrixInv](@ref pm_matrixInv).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_matrixInv) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
    use pm_matrixCopy, only: getMatCopy, lfpack, dia
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_matrixLUP, only: setMatLUP
   !use pm_matrixInit, only: setMatInit
    use pm_matrixCopy, only: setMatCopy, rdpack
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatInv_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IMP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInvDef_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInvDef_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInvDef_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInvDef_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInvDef_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInvDef_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInvDef_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInvDef_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInvDef_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInvDef_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Det_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInvDet_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInvDet_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInvDet_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInvDet_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInvDet_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInvDet_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInvDet_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInvDet_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInvDet_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInvDet_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Det_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Inf_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInvInf_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInvInf_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInvInf_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInvInf_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInvInf_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInvInf_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInvInf_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInvInf_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInvInf_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInvInf_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Inf_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInvCUD_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInvCUD_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInvCUD_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInvCUD_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInvCUD_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInvCUD_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInvCUD_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInvCUD_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInvCUD_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInvCUD_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInvCLD_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInvCLD_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInvCLD_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInvCLD_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInvCLD_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInvCLD_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInvCLD_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInvCLD_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInvCLD_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInvCLD_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInvCUU_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInvCUU_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInvCUU_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInvCUU_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInvCUU_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInvCUU_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInvCUU_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInvCUU_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInvCUU_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInvCUU_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInvCLU_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInvCLU_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInvCLU_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInvCLU_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInvCLU_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInvCLU_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInvCLU_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInvCLU_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInvCLU_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInvCLU_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LUP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInvLUP_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInvLUP_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInvLUP_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInvLUP_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInvLUP_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInvLUP_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInvLUP_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInvLUP_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInvLUP_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInvLUP_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LUP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CCU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInvCCU_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInvCCU_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInvCCU_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInvCCU_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInvCCU_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInvCCU_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInvCCU_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInvCCU_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInvCCU_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInvCCU_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CCU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CCL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatInvCCL_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatInvCCL_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatInvCCL_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatInvCCL_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatInvCCL_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatInvCCL_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatInvCCL_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatInvCCL_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatInvCCL_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatInvCCL_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CCL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IMP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMatInv_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatInv_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IMP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FUL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInvCUD_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInvCUD_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInvCUD_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInvCUD_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInvCUD_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInvCUD_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInvCUD_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInvCUD_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInvCUD_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInvCUD_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInvCLD_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInvCLD_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInvCLD_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInvCLD_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInvCLD_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInvCLD_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInvCLD_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInvCLD_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInvCLD_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInvCLD_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInvCUU_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInvCUU_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInvCUU_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInvCUU_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInvCUU_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInvCUU_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInvCUU_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInvCUU_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInvCUU_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInvCUU_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInvCLU_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInvCLU_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInvCLU_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInvCLU_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInvCLU_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInvCLU_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInvCLU_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInvCLU_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInvCLU_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInvCLU_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LUP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInvLUP_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInvLUP_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInvLUP_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInvLUP_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInvLUP_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInvLUP_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInvLUP_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInvLUP_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInvLUP_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInvLUP_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LUP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CCU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInvCCU_FUL_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInvCCU_FUL_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInvCCU_FUL_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInvCCU_FUL_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInvCCU_FUL_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInvCCU_FUL_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInvCCU_FUL_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInvCCU_FUL_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInvCCU_FUL_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInvCCU_FUL_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CCU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CCL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInvCCL_FUL_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInvCCL_FUL_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInvCCL_FUL_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInvCCL_FUL_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInvCCL_FUL_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInvCCL_FUL_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInvCCL_FUL_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInvCCL_FUL_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInvCCL_FUL_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInvCCL_FUL_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CCL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FUL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IMP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatInv_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatInv_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IMP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CCU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInvCCU_UXD_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInvCCU_UXD_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInvCCU_UXD_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInvCCU_UXD_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInvCCU_UXD_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInvCCU_UXD_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInvCCU_UXD_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInvCCU_UXD_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInvCCU_UXD_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInvCCU_UXD_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CCU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CCL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInvCCL_UXD_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInvCCL_UXD_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInvCCL_UXD_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInvCCL_UXD_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInvCCL_UXD_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInvCCL_UXD_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInvCCL_UXD_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInvCCL_UXD_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInvCCL_UXD_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInvCCL_UXD_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CCL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IMP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatInv_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatInv_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IMP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CCU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInvCCU_XLD_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInvCCU_XLD_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInvCCU_XLD_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInvCCU_XLD_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInvCCU_XLD_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInvCCU_XLD_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInvCCU_XLD_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInvCCU_XLD_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInvCCU_XLD_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInvCCU_XLD_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CCU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CCL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatInvCCL_XLD_IMP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatInvCCL_XLD_IMP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatInvCCL_XLD_IMP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatInvCCL_XLD_IMP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatInvCCL_XLD_IMP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatInvCCL_XLD_IMP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatInvCCL_XLD_IMP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatInvCCL_XLD_IMP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixInv@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatInvCCL_XLD_IMP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatInvCCL_XLD_IMP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixInv@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CCL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IMP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatInv_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines