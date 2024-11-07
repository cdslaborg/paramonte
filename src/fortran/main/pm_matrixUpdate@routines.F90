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
!>  This file contains procedure implementations of [pm_matrixUpdate](@ref pm_matrixUpdate).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_matrixUpdate) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    ! \bug Bypass Intel `ifort` 2022 compiler bug for too many use statements in submodule procedures.
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif
    use pm_blas, only: blasSYRK, blasHERK

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatUpdateR1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Fixed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatUpdateR1F_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatUpdateR1F_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatUpdateR1F_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatUpdateR1F_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatUpdateR1F_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CHM_ENABLED 1
#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatUpdateR1H_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatUpdateR1H_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatUpdateR1H_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatUpdateR1H_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatUpdateR1H_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED
#undef CHM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1
#define CSM_ENABLED 1

#if CK5_ENABLED
    module procedure setMatUpdateR1F_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatUpdateR1F_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatUpdateR1F_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatUpdateR1F_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatUpdateR1F_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef CSM_ENABLED
#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatUpdateR1F_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatUpdateR1F_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatUpdateR1F_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatUpdateR1F_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatUpdateR1F_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Fixed_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Alpha_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatUpdateR1A_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatUpdateR1A_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatUpdateR1A_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatUpdateR1A_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatUpdateR1A_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1
#define CHM_ENABLED 1

#if CK5_ENABLED
    module procedure setMatUpdateR1AH_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatUpdateR1AH_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatUpdateR1AH_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatUpdateR1AH_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatUpdateR1AH_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef CHM_ENABLED
#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1
#define CSM_ENABLED 1

#if CK5_ENABLED
    module procedure setMatUpdateR1A_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatUpdateR1A_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatUpdateR1A_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatUpdateR1A_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatUpdateR1A_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef CSM_ENABLED
#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatUpdateR1A_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatUpdateR1A_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatUpdateR1A_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatUpdateR1A_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatUpdateR1A_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Alpha_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatUpdateR1_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatUpdate_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define shrk_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ASS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CSM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure shrk_ASS_CSM_SLD_ONO_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure shrk_ASS_CSM_SLD_ONO_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure shrk_ASS_CSM_SLD_ONO_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure shrk_ASS_CSM_SLD_ONO_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure shrk_ASS_CSM_SLD_ONO_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure shrk_ASS_CSM_SLD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure shrk_ASS_CSM_SLD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure shrk_ASS_CSM_SLD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure shrk_ASS_CSM_SLD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure shrk_ASS_CSM_SLD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure shrk_ASS_CSM_SLD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure shrk_ASS_CSM_SLD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure shrk_ASS_CSM_SLD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure shrk_ASS_CSM_SLD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure shrk_ASS_CSM_SLD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure shrk_ASS_CSM_SLD_OTP_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure shrk_ASS_CSM_SLD_OTP_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure shrk_ASS_CSM_SLD_OTP_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure shrk_ASS_CSM_SLD_OTP_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure shrk_ASS_CSM_SLD_OTP_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure shrk_ASS_CSM_SLD_OTP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure shrk_ASS_CSM_SLD_OTP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure shrk_ASS_CSM_SLD_OTP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure shrk_ASS_CSM_SLD_OTP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure shrk_ASS_CSM_SLD_OTP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure shrk_ASS_CSM_SLD_OTP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure shrk_ASS_CSM_SLD_OTP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure shrk_ASS_CSM_SLD_OTP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure shrk_ASS_CSM_SLD_OTP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure shrk_ASS_CSM_SLD_OTP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure shrk_ASS_CSM_SUD_ONO_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure shrk_ASS_CSM_SUD_ONO_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure shrk_ASS_CSM_SUD_ONO_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure shrk_ASS_CSM_SUD_ONO_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure shrk_ASS_CSM_SUD_ONO_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure shrk_ASS_CSM_SUD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure shrk_ASS_CSM_SUD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure shrk_ASS_CSM_SUD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure shrk_ASS_CSM_SUD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure shrk_ASS_CSM_SUD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure shrk_ASS_CSM_SUD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure shrk_ASS_CSM_SUD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure shrk_ASS_CSM_SUD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure shrk_ASS_CSM_SUD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure shrk_ASS_CSM_SUD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure shrk_ASS_CSM_SUD_OTP_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure shrk_ASS_CSM_SUD_OTP_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure shrk_ASS_CSM_SUD_OTP_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure shrk_ASS_CSM_SUD_OTP_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure shrk_ASS_CSM_SUD_OTP_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure shrk_ASS_CSM_SUD_OTP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure shrk_ASS_CSM_SUD_OTP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure shrk_ASS_CSM_SUD_OTP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure shrk_ASS_CSM_SUD_OTP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure shrk_ASS_CSM_SUD_OTP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure shrk_ASS_CSM_SUD_OTP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure shrk_ASS_CSM_SUD_OTP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure shrk_ASS_CSM_SUD_OTP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure shrk_ASS_CSM_SUD_OTP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure shrk_ASS_CSM_SUD_OTP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CSM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CHM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure shrk_ASS_CHM_SLD_ONO_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure shrk_ASS_CHM_SLD_ONO_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure shrk_ASS_CHM_SLD_ONO_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure shrk_ASS_CHM_SLD_ONO_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure shrk_ASS_CHM_SLD_ONO_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure shrk_ASS_CHM_SLD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure shrk_ASS_CHM_SLD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure shrk_ASS_CHM_SLD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure shrk_ASS_CHM_SLD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure shrk_ASS_CHM_SLD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure shrk_ASS_CHM_SLD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure shrk_ASS_CHM_SLD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure shrk_ASS_CHM_SLD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure shrk_ASS_CHM_SLD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure shrk_ASS_CHM_SLD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure shrk_ASS_CHM_SLD_OTP_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure shrk_ASS_CHM_SLD_OTP_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure shrk_ASS_CHM_SLD_OTP_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure shrk_ASS_CHM_SLD_OTP_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure shrk_ASS_CHM_SLD_OTP_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure shrk_ASS_CHM_SLD_OTP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure shrk_ASS_CHM_SLD_OTP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure shrk_ASS_CHM_SLD_OTP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure shrk_ASS_CHM_SLD_OTP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure shrk_ASS_CHM_SLD_OTP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure shrk_ASS_CHM_SLD_OTP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure shrk_ASS_CHM_SLD_OTP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure shrk_ASS_CHM_SLD_OTP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure shrk_ASS_CHM_SLD_OTP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure shrk_ASS_CHM_SLD_OTP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure shrk_ASS_CHM_SUD_ONO_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure shrk_ASS_CHM_SUD_ONO_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure shrk_ASS_CHM_SUD_ONO_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure shrk_ASS_CHM_SUD_ONO_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure shrk_ASS_CHM_SUD_ONO_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure shrk_ASS_CHM_SUD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure shrk_ASS_CHM_SUD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure shrk_ASS_CHM_SUD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure shrk_ASS_CHM_SUD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure shrk_ASS_CHM_SUD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure shrk_ASS_CHM_SUD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure shrk_ASS_CHM_SUD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure shrk_ASS_CHM_SUD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure shrk_ASS_CHM_SUD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure shrk_ASS_CHM_SUD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure shrk_ASS_CHM_SUD_OTP_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure shrk_ASS_CHM_SUD_OTP_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure shrk_ASS_CHM_SUD_OTP_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure shrk_ASS_CHM_SUD_OTP_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure shrk_ASS_CHM_SUD_OTP_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure shrk_ASS_CHM_SUD_OTP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure shrk_ASS_CHM_SUD_OTP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure shrk_ASS_CHM_SUD_OTP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure shrk_ASS_CHM_SUD_OTP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure shrk_ASS_CHM_SUD_OTP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure shrk_ASS_CHM_SUD_OTP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure shrk_ASS_CHM_SUD_OTP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure shrk_ASS_CHM_SUD_OTP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure shrk_ASS_CHM_SUD_OTP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure shrk_ASS_CHM_SUD_OTP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ASS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define EXP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CSM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure shrk_EXP_CSM_SLD_ONO_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure shrk_EXP_CSM_SLD_ONO_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure shrk_EXP_CSM_SLD_ONO_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure shrk_EXP_CSM_SLD_ONO_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure shrk_EXP_CSM_SLD_ONO_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure shrk_EXP_CSM_SLD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure shrk_EXP_CSM_SLD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure shrk_EXP_CSM_SLD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure shrk_EXP_CSM_SLD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure shrk_EXP_CSM_SLD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure shrk_EXP_CSM_SLD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure shrk_EXP_CSM_SLD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure shrk_EXP_CSM_SLD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure shrk_EXP_CSM_SLD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure shrk_EXP_CSM_SLD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure shrk_EXP_CSM_SLD_OTP_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure shrk_EXP_CSM_SLD_OTP_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure shrk_EXP_CSM_SLD_OTP_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure shrk_EXP_CSM_SLD_OTP_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure shrk_EXP_CSM_SLD_OTP_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure shrk_EXP_CSM_SLD_OTP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure shrk_EXP_CSM_SLD_OTP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure shrk_EXP_CSM_SLD_OTP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure shrk_EXP_CSM_SLD_OTP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure shrk_EXP_CSM_SLD_OTP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure shrk_EXP_CSM_SLD_OTP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure shrk_EXP_CSM_SLD_OTP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure shrk_EXP_CSM_SLD_OTP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure shrk_EXP_CSM_SLD_OTP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure shrk_EXP_CSM_SLD_OTP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure shrk_EXP_CSM_SUD_ONO_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure shrk_EXP_CSM_SUD_ONO_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure shrk_EXP_CSM_SUD_ONO_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure shrk_EXP_CSM_SUD_ONO_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure shrk_EXP_CSM_SUD_ONO_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure shrk_EXP_CSM_SUD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure shrk_EXP_CSM_SUD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure shrk_EXP_CSM_SUD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure shrk_EXP_CSM_SUD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure shrk_EXP_CSM_SUD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure shrk_EXP_CSM_SUD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure shrk_EXP_CSM_SUD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure shrk_EXP_CSM_SUD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure shrk_EXP_CSM_SUD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure shrk_EXP_CSM_SUD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure shrk_EXP_CSM_SUD_OTP_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure shrk_EXP_CSM_SUD_OTP_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure shrk_EXP_CSM_SUD_OTP_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure shrk_EXP_CSM_SUD_OTP_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure shrk_EXP_CSM_SUD_OTP_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure shrk_EXP_CSM_SUD_OTP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure shrk_EXP_CSM_SUD_OTP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure shrk_EXP_CSM_SUD_OTP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure shrk_EXP_CSM_SUD_OTP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure shrk_EXP_CSM_SUD_OTP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure shrk_EXP_CSM_SUD_OTP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure shrk_EXP_CSM_SUD_OTP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure shrk_EXP_CSM_SUD_OTP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure shrk_EXP_CSM_SUD_OTP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure shrk_EXP_CSM_SUD_OTP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CSM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CHM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure shrk_EXP_CHM_SLD_ONO_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure shrk_EXP_CHM_SLD_ONO_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure shrk_EXP_CHM_SLD_ONO_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure shrk_EXP_CHM_SLD_ONO_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure shrk_EXP_CHM_SLD_ONO_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure shrk_EXP_CHM_SLD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure shrk_EXP_CHM_SLD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure shrk_EXP_CHM_SLD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure shrk_EXP_CHM_SLD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure shrk_EXP_CHM_SLD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure shrk_EXP_CHM_SLD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure shrk_EXP_CHM_SLD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure shrk_EXP_CHM_SLD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure shrk_EXP_CHM_SLD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure shrk_EXP_CHM_SLD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure shrk_EXP_CHM_SLD_OTP_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure shrk_EXP_CHM_SLD_OTP_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure shrk_EXP_CHM_SLD_OTP_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure shrk_EXP_CHM_SLD_OTP_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure shrk_EXP_CHM_SLD_OTP_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure shrk_EXP_CHM_SLD_OTP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure shrk_EXP_CHM_SLD_OTP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure shrk_EXP_CHM_SLD_OTP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure shrk_EXP_CHM_SLD_OTP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure shrk_EXP_CHM_SLD_OTP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure shrk_EXP_CHM_SLD_OTP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure shrk_EXP_CHM_SLD_OTP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure shrk_EXP_CHM_SLD_OTP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure shrk_EXP_CHM_SLD_OTP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure shrk_EXP_CHM_SLD_OTP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SUD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure shrk_EXP_CHM_SUD_ONO_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure shrk_EXP_CHM_SUD_ONO_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure shrk_EXP_CHM_SUD_ONO_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure shrk_EXP_CHM_SUD_ONO_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure shrk_EXP_CHM_SUD_ONO_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure shrk_EXP_CHM_SUD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure shrk_EXP_CHM_SUD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure shrk_EXP_CHM_SUD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure shrk_EXP_CHM_SUD_ONO_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure shrk_EXP_CHM_SUD_ONO_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure shrk_EXP_CHM_SUD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure shrk_EXP_CHM_SUD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure shrk_EXP_CHM_SUD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure shrk_EXP_CHM_SUD_ONO_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure shrk_EXP_CHM_SUD_ONO_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure shrk_EXP_CHM_SUD_OTP_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure shrk_EXP_CHM_SUD_OTP_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure shrk_EXP_CHM_SUD_OTP_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure shrk_EXP_CHM_SUD_OTP_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure shrk_EXP_CHM_SUD_OTP_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure shrk_EXP_CHM_SUD_OTP_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure shrk_EXP_CHM_SUD_OTP_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure shrk_EXP_CHM_SUD_OTP_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure shrk_EXP_CHM_SUD_OTP_CK2
        use pm_kind, only: TKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure shrk_EXP_CHM_SUD_OTP_CK1
        use pm_kind, only: TKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure shrk_EXP_CHM_SUD_OTP_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure shrk_EXP_CHM_SUD_OTP_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure shrk_EXP_CHM_SUD_OTP_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure shrk_EXP_CHM_SUD_OTP_RK2
        use pm_kind, only: TKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure shrk_EXP_CHM_SUD_OTP_RK1
        use pm_kind, only: TKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixUpdate@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SUD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef EXP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef shrk_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatUpdate_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatUpdateTriang_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define lowDiaC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatUpdateTriangCSOLAS_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatUpdateTriangCSOLAS_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatUpdateTriangCSOLAS_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatUpdateTriangCSOLAS_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatUpdateTriangCSOLAS_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatUpdateTriangCSOLAS_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatUpdateTriangCSOLAS_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatUpdateTriangCSOLAS_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatUpdateTriangCSOLAS_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatUpdateTriangCSOLAS_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatUpdateTriangCSOLAS_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatUpdateTriangCSOLAS_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatUpdateTriangCSOLAS_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatUpdateTriangCSOLAS_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatUpdateTriangCSOLAS_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatUpdateTriangCSOLSA_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatUpdateTriangCSOLSA_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatUpdateTriangCSOLSA_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatUpdateTriangCSOLSA_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatUpdateTriangCSOLSA_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatUpdateTriangCSOLSA_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatUpdateTriangCSOLSA_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatUpdateTriangCSOLSA_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatUpdateTriangCSOLSA_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatUpdateTriangCSOLSA_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatUpdateTriangCSOLSA_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatUpdateTriangCSOLSA_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatUpdateTriangCSOLSA_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatUpdateTriangCSOLSA_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatUpdateTriangCSOLSA_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AH_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatUpdateTriangCSOLAH_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatUpdateTriangCSOLAH_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatUpdateTriangCSOLAH_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatUpdateTriangCSOLAH_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatUpdateTriangCSOLAH_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AH_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define HA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatUpdateTriangCSOLHA_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatUpdateTriangCSOLHA_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatUpdateTriangCSOLHA_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatUpdateTriangCSOLHA_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatUpdateTriangCSOLHA_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef HA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef lowDiaC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define uppDiaC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatUpdateTriangCSOUAS_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatUpdateTriangCSOUAS_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatUpdateTriangCSOUAS_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatUpdateTriangCSOUAS_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatUpdateTriangCSOUAS_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatUpdateTriangCSOUAS_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatUpdateTriangCSOUAS_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatUpdateTriangCSOUAS_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatUpdateTriangCSOUAS_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatUpdateTriangCSOUAS_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatUpdateTriangCSOUAS_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatUpdateTriangCSOUAS_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatUpdateTriangCSOUAS_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatUpdateTriangCSOUAS_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatUpdateTriangCSOUAS_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMatUpdateTriangCSOUSA_IK5
        use pm_kind, only: TKG => IK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMatUpdateTriangCSOUSA_IK4
        use pm_kind, only: TKG => IK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMatUpdateTriangCSOUSA_IK3
        use pm_kind, only: TKG => IK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMatUpdateTriangCSOUSA_IK2
        use pm_kind, only: TKG => IK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMatUpdateTriangCSOUSA_IK1
        use pm_kind, only: TKG => IK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatUpdateTriangCSOUSA_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatUpdateTriangCSOUSA_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatUpdateTriangCSOUSA_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatUpdateTriangCSOUSA_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatUpdateTriangCSOUSA_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatUpdateTriangCSOUSA_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatUpdateTriangCSOUSA_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatUpdateTriangCSOUSA_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatUpdateTriangCSOUSA_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatUpdateTriangCSOUSA_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AH_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatUpdateTriangCSOUAH_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatUpdateTriangCSOUAH_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatUpdateTriangCSOUAH_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatUpdateTriangCSOUAH_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatUpdateTriangCSOUAH_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AH_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define HA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatUpdateTriangCSOUHA_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatUpdateTriangCSOUHA_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatUpdateTriangCSOUHA_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatUpdateTriangCSOUHA_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatUpdateTriangCSOUHA_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixUpdate@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef HA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef uppDiaC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatUpdateTriang_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines