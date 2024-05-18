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
!>  This file contains procedure implementations of [pm_container](@ref pm_container).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Monday 2:21 AM, August 30, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_container) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_logicalCompare, only: operator(<), operator(>), operator(<=), operator(>=), operator(==), operator(/=)
    use pm_complexCompareLex, only: operator(<), operator(>), operator(<=), operator(>=)
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define constructCon_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1
    module procedure css_typer_D0
        use pm_kind, only: SKC => SK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  SK_ENABLED

#define IK_ENABLED 1
    module procedure csi_typer_D0
        use pm_kind, only: IKC => IK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  IK_ENABLED

#define LK_ENABLED 1
    module procedure csl_typer_D0
        use pm_kind, only: LKC => LK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  LK_ENABLED

#define CK_ENABLED 1
    module procedure csc_typer_D0
        use pm_kind, only: CKC => CK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  CK_ENABLED

#define RK_ENABLED 1
    module procedure csr_typer_D0
        use pm_kind, only: RKC => RK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  RK_ENABLED

#define PK_ENABLED 1
    module procedure csp_typer_D0
#include "pm_container@routines.inc.F90"
    end procedure
#undef  PK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef BS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure constructCon_D0_PSSK5
        use pm_kind, only: SKC => SK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure constructCon_D0_PSSK4
        use pm_kind, only: SKC => SK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure constructCon_D0_PSSK3
        use pm_kind, only: SKC => SK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure constructCon_D0_PSSK2
        use pm_kind, only: SKC => SK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure constructCon_D0_PSSK1
        use pm_kind, only: SKC => SK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure constructCon_D0_PSIK5
        use pm_kind, only: IKC => IK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure constructCon_D0_PSIK4
        use pm_kind, only: IKC => IK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure constructCon_D0_PSIK3
        use pm_kind, only: IKC => IK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure constructCon_D0_PSIK2
        use pm_kind, only: IKC => IK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure constructCon_D0_PSIK1
        use pm_kind, only: IKC => IK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure constructCon_D0_PSLK5
        use pm_kind, only: LKC => LK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure constructCon_D0_PSLK4
        use pm_kind, only: LKC => LK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure constructCon_D0_PSLK3
        use pm_kind, only: LKC => LK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure constructCon_D0_PSLK2
        use pm_kind, only: LKC => LK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure constructCon_D0_PSLK1
        use pm_kind, only: LKC => LK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure constructCon_D0_PSCK5
        use pm_kind, only: CKC => CK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure constructCon_D0_PSCK4
        use pm_kind, only: CKC => CK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure constructCon_D0_PSCK3
        use pm_kind, only: CKC => CK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure constructCon_D0_PSCK2
        use pm_kind, only: CKC => CK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure constructCon_D0_PSCK1
        use pm_kind, only: CKC => CK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure constructCon_D0_PSRK5
        use pm_kind, only: RKC => RK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure constructCon_D0_PSRK4
        use pm_kind, only: RKC => RK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure constructCon_D0_PSRK3
        use pm_kind, only: RKC => RK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure constructCon_D0_PSRK2
        use pm_kind, only: RKC => RK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure constructCon_D0_PSRK1
        use pm_kind, only: RKC => RK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PS_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef constructCon_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isless_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1
    module procedure isless_D0_D0_BSSK
        use pm_kind, only: SKC => SK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  SK_ENABLED

#define IK_ENABLED 1
    module procedure isless_D0_D0_BSIK
        use pm_kind, only: IKC => IK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  IK_ENABLED

#define LK_ENABLED 1
    module procedure isless_D0_D0_BSLK
        use pm_kind, only: LKC => LK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  LK_ENABLED

#define CK_ENABLED 1
    module procedure isless_D0_D0_BSCK
        use pm_kind, only: CKC => CK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  CK_ENABLED

#define RK_ENABLED 1
    module procedure isless_D0_D0_BSRK
        use pm_kind, only: RKC => RK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef BS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isless_D0_D0_PSSK5
        use pm_kind, only: SKC => SK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isless_D0_D0_PSSK4
        use pm_kind, only: SKC => SK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isless_D0_D0_PSSK3
        use pm_kind, only: SKC => SK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isless_D0_D0_PSSK2
        use pm_kind, only: SKC => SK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isless_D0_D0_PSSK1
        use pm_kind, only: SKC => SK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure isless_D0_D0_PSIK5
        use pm_kind, only: IKC => IK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure isless_D0_D0_PSIK4
        use pm_kind, only: IKC => IK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure isless_D0_D0_PSIK3
        use pm_kind, only: IKC => IK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure isless_D0_D0_PSIK2
        use pm_kind, only: IKC => IK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure isless_D0_D0_PSIK1
        use pm_kind, only: IKC => IK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure isless_D0_D0_PSLK5
        use pm_kind, only: LKC => LK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure isless_D0_D0_PSLK4
        use pm_kind, only: LKC => LK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure isless_D0_D0_PSLK3
        use pm_kind, only: LKC => LK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure isless_D0_D0_PSLK2
        use pm_kind, only: LKC => LK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure isless_D0_D0_PSLK1
        use pm_kind, only: LKC => LK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isless_D0_D0_PSCK5
        use pm_kind, only: CKC => CK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isless_D0_D0_PSCK4
        use pm_kind, only: CKC => CK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isless_D0_D0_PSCK3
        use pm_kind, only: CKC => CK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isless_D0_D0_PSCK2
        use pm_kind, only: CKC => CK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isless_D0_D0_PSCK1
        use pm_kind, only: CKC => CK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isless_D0_D0_PSRK5
        use pm_kind, only: RKC => RK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isless_D0_D0_PSRK4
        use pm_kind, only: RKC => RK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isless_D0_D0_PSRK3
        use pm_kind, only: RKC => RK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isless_D0_D0_PSRK2
        use pm_kind, only: RKC => RK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isless_D0_D0_PSRK1
        use pm_kind, only: RKC => RK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PS_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isless_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ismore_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1
    module procedure ismore_D0_D0_BSSK
        use pm_kind, only: SKC => SK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  SK_ENABLED

#define IK_ENABLED 1
    module procedure ismore_D0_D0_BSIK
        use pm_kind, only: IKC => IK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  IK_ENABLED

#define LK_ENABLED 1
    module procedure ismore_D0_D0_BSLK
        use pm_kind, only: LKC => LK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  LK_ENABLED

#define CK_ENABLED 1
    module procedure ismore_D0_D0_BSCK
        use pm_kind, only: CKC => CK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  CK_ENABLED

#define RK_ENABLED 1
    module procedure ismore_D0_D0_BSRK
        use pm_kind, only: RKC => RK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef BS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure ismore_D0_D0_PSSK5
        use pm_kind, only: SKC => SK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure ismore_D0_D0_PSSK4
        use pm_kind, only: SKC => SK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure ismore_D0_D0_PSSK3
        use pm_kind, only: SKC => SK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure ismore_D0_D0_PSSK2
        use pm_kind, only: SKC => SK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure ismore_D0_D0_PSSK1
        use pm_kind, only: SKC => SK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure ismore_D0_D0_PSIK5
        use pm_kind, only: IKC => IK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure ismore_D0_D0_PSIK4
        use pm_kind, only: IKC => IK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure ismore_D0_D0_PSIK3
        use pm_kind, only: IKC => IK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure ismore_D0_D0_PSIK2
        use pm_kind, only: IKC => IK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure ismore_D0_D0_PSIK1
        use pm_kind, only: IKC => IK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure ismore_D0_D0_PSLK5
        use pm_kind, only: LKC => LK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure ismore_D0_D0_PSLK4
        use pm_kind, only: LKC => LK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure ismore_D0_D0_PSLK3
        use pm_kind, only: LKC => LK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure ismore_D0_D0_PSLK2
        use pm_kind, only: LKC => LK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure ismore_D0_D0_PSLK1
        use pm_kind, only: LKC => LK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure ismore_D0_D0_PSCK5
        use pm_kind, only: CKC => CK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure ismore_D0_D0_PSCK4
        use pm_kind, only: CKC => CK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure ismore_D0_D0_PSCK3
        use pm_kind, only: CKC => CK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure ismore_D0_D0_PSCK2
        use pm_kind, only: CKC => CK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure ismore_D0_D0_PSCK1
        use pm_kind, only: CKC => CK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure ismore_D0_D0_PSRK5
        use pm_kind, only: RKC => RK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure ismore_D0_D0_PSRK4
        use pm_kind, only: RKC => RK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure ismore_D0_D0_PSRK3
        use pm_kind, only: RKC => RK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure ismore_D0_D0_PSRK2
        use pm_kind, only: RKC => RK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure ismore_D0_D0_PSRK1
        use pm_kind, only: RKC => RK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PS_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ismore_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isleq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1
    module procedure isleq_D0_D0_BSSK
        use pm_kind, only: SKC => SK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  SK_ENABLED

#define IK_ENABLED 1
    module procedure isleq_D0_D0_BSIK
        use pm_kind, only: IKC => IK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  IK_ENABLED

#define LK_ENABLED 1
    module procedure isleq_D0_D0_BSLK
        use pm_kind, only: LKC => LK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  LK_ENABLED

#define CK_ENABLED 1
    module procedure isleq_D0_D0_BSCK
        use pm_kind, only: CKC => CK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  CK_ENABLED

#define RK_ENABLED 1
    module procedure isleq_D0_D0_BSRK
        use pm_kind, only: RKC => RK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef BS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isleq_D0_D0_PSSK5
        use pm_kind, only: SKC => SK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isleq_D0_D0_PSSK4
        use pm_kind, only: SKC => SK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isleq_D0_D0_PSSK3
        use pm_kind, only: SKC => SK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isleq_D0_D0_PSSK2
        use pm_kind, only: SKC => SK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isleq_D0_D0_PSSK1
        use pm_kind, only: SKC => SK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure isleq_D0_D0_PSIK5
        use pm_kind, only: IKC => IK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure isleq_D0_D0_PSIK4
        use pm_kind, only: IKC => IK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure isleq_D0_D0_PSIK3
        use pm_kind, only: IKC => IK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure isleq_D0_D0_PSIK2
        use pm_kind, only: IKC => IK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure isleq_D0_D0_PSIK1
        use pm_kind, only: IKC => IK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure isleq_D0_D0_PSLK5
        use pm_kind, only: LKC => LK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure isleq_D0_D0_PSLK4
        use pm_kind, only: LKC => LK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure isleq_D0_D0_PSLK3
        use pm_kind, only: LKC => LK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure isleq_D0_D0_PSLK2
        use pm_kind, only: LKC => LK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure isleq_D0_D0_PSLK1
        use pm_kind, only: LKC => LK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isleq_D0_D0_PSCK5
        use pm_kind, only: CKC => CK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isleq_D0_D0_PSCK4
        use pm_kind, only: CKC => CK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isleq_D0_D0_PSCK3
        use pm_kind, only: CKC => CK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isleq_D0_D0_PSCK2
        use pm_kind, only: CKC => CK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isleq_D0_D0_PSCK1
        use pm_kind, only: CKC => CK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isleq_D0_D0_PSRK5
        use pm_kind, only: RKC => RK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isleq_D0_D0_PSRK4
        use pm_kind, only: RKC => RK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isleq_D0_D0_PSRK3
        use pm_kind, only: RKC => RK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isleq_D0_D0_PSRK2
        use pm_kind, only: RKC => RK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isleq_D0_D0_PSRK1
        use pm_kind, only: RKC => RK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PS_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isleq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ismeq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1
    module procedure ismeq_D0_D0_BSSK
        use pm_kind, only: SKC => SK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  SK_ENABLED

#define IK_ENABLED 1
    module procedure ismeq_D0_D0_BSIK
        use pm_kind, only: IKC => IK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  IK_ENABLED

#define LK_ENABLED 1
    module procedure ismeq_D0_D0_BSLK
        use pm_kind, only: LKC => LK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  LK_ENABLED

#define CK_ENABLED 1
    module procedure ismeq_D0_D0_BSCK
        use pm_kind, only: CKC => CK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  CK_ENABLED

#define RK_ENABLED 1
    module procedure ismeq_D0_D0_BSRK
        use pm_kind, only: RKC => RK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef BS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure ismeq_D0_D0_PSSK5
        use pm_kind, only: SKC => SK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure ismeq_D0_D0_PSSK4
        use pm_kind, only: SKC => SK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure ismeq_D0_D0_PSSK3
        use pm_kind, only: SKC => SK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure ismeq_D0_D0_PSSK2
        use pm_kind, only: SKC => SK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure ismeq_D0_D0_PSSK1
        use pm_kind, only: SKC => SK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure ismeq_D0_D0_PSIK5
        use pm_kind, only: IKC => IK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure ismeq_D0_D0_PSIK4
        use pm_kind, only: IKC => IK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure ismeq_D0_D0_PSIK3
        use pm_kind, only: IKC => IK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure ismeq_D0_D0_PSIK2
        use pm_kind, only: IKC => IK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure ismeq_D0_D0_PSIK1
        use pm_kind, only: IKC => IK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure ismeq_D0_D0_PSLK5
        use pm_kind, only: LKC => LK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure ismeq_D0_D0_PSLK4
        use pm_kind, only: LKC => LK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure ismeq_D0_D0_PSLK3
        use pm_kind, only: LKC => LK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure ismeq_D0_D0_PSLK2
        use pm_kind, only: LKC => LK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure ismeq_D0_D0_PSLK1
        use pm_kind, only: LKC => LK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure ismeq_D0_D0_PSCK5
        use pm_kind, only: CKC => CK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure ismeq_D0_D0_PSCK4
        use pm_kind, only: CKC => CK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure ismeq_D0_D0_PSCK3
        use pm_kind, only: CKC => CK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure ismeq_D0_D0_PSCK2
        use pm_kind, only: CKC => CK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure ismeq_D0_D0_PSCK1
        use pm_kind, only: CKC => CK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure ismeq_D0_D0_PSRK5
        use pm_kind, only: RKC => RK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure ismeq_D0_D0_PSRK4
        use pm_kind, only: RKC => RK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure ismeq_D0_D0_PSRK3
        use pm_kind, only: RKC => RK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure ismeq_D0_D0_PSRK2
        use pm_kind, only: RKC => RK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure ismeq_D0_D0_PSRK1
        use pm_kind, only: RKC => RK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PS_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ismeq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isneq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1
    module procedure isneq_D0_D0_BSSK
        use pm_kind, only: SKC => SK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  SK_ENABLED

#define IK_ENABLED 1
    module procedure isneq_D0_D0_BSIK
        use pm_kind, only: IKC => IK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  IK_ENABLED

#define LK_ENABLED 1
    module procedure isneq_D0_D0_BSLK
        use pm_kind, only: LKC => LK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  LK_ENABLED

#define CK_ENABLED 1
    module procedure isneq_D0_D0_BSCK
        use pm_kind, only: CKC => CK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  CK_ENABLED

#define RK_ENABLED 1
    module procedure isneq_D0_D0_BSRK
        use pm_kind, only: RKC => RK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef BS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isneq_D0_D0_PSSK5
        use pm_kind, only: SKC => SK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isneq_D0_D0_PSSK4
        use pm_kind, only: SKC => SK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isneq_D0_D0_PSSK3
        use pm_kind, only: SKC => SK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isneq_D0_D0_PSSK2
        use pm_kind, only: SKC => SK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isneq_D0_D0_PSSK1
        use pm_kind, only: SKC => SK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure isneq_D0_D0_PSIK5
        use pm_kind, only: IKC => IK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure isneq_D0_D0_PSIK4
        use pm_kind, only: IKC => IK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure isneq_D0_D0_PSIK3
        use pm_kind, only: IKC => IK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure isneq_D0_D0_PSIK2
        use pm_kind, only: IKC => IK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure isneq_D0_D0_PSIK1
        use pm_kind, only: IKC => IK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure isneq_D0_D0_PSLK5
        use pm_kind, only: LKC => LK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure isneq_D0_D0_PSLK4
        use pm_kind, only: LKC => LK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure isneq_D0_D0_PSLK3
        use pm_kind, only: LKC => LK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure isneq_D0_D0_PSLK2
        use pm_kind, only: LKC => LK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure isneq_D0_D0_PSLK1
        use pm_kind, only: LKC => LK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure isneq_D0_D0_PSCK5
        use pm_kind, only: CKC => CK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isneq_D0_D0_PSCK4
        use pm_kind, only: CKC => CK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isneq_D0_D0_PSCK3
        use pm_kind, only: CKC => CK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isneq_D0_D0_PSCK2
        use pm_kind, only: CKC => CK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isneq_D0_D0_PSCK1
        use pm_kind, only: CKC => CK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isneq_D0_D0_PSRK5
        use pm_kind, only: RKC => RK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isneq_D0_D0_PSRK4
        use pm_kind, only: RKC => RK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isneq_D0_D0_PSRK3
        use pm_kind, only: RKC => RK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isneq_D0_D0_PSRK2
        use pm_kind, only: RKC => RK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isneq_D0_D0_PSRK1
        use pm_kind, only: RKC => RK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PS_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isneq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define iseq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1
    module procedure iseq_D0_D0_BSSK
        use pm_kind, only: SKC => SK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  SK_ENABLED

#define IK_ENABLED 1
    module procedure iseq_D0_D0_BSIK
        use pm_kind, only: IKC => IK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  IK_ENABLED

#define LK_ENABLED 1
    module procedure iseq_D0_D0_BSLK
        use pm_kind, only: LKC => LK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  LK_ENABLED

#define CK_ENABLED 1
    module procedure iseq_D0_D0_BSCK
        use pm_kind, only: CKC => CK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  CK_ENABLED

#define RK_ENABLED 1
    module procedure iseq_D0_D0_BSRK
        use pm_kind, only: RKC => RK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef BS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure iseq_D0_D0_PSSK5
        use pm_kind, only: SKC => SK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure iseq_D0_D0_PSSK4
        use pm_kind, only: SKC => SK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure iseq_D0_D0_PSSK3
        use pm_kind, only: SKC => SK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure iseq_D0_D0_PSSK2
        use pm_kind, only: SKC => SK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure iseq_D0_D0_PSSK1
        use pm_kind, only: SKC => SK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure iseq_D0_D0_PSIK5
        use pm_kind, only: IKC => IK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure iseq_D0_D0_PSIK4
        use pm_kind, only: IKC => IK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure iseq_D0_D0_PSIK3
        use pm_kind, only: IKC => IK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure iseq_D0_D0_PSIK2
        use pm_kind, only: IKC => IK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure iseq_D0_D0_PSIK1
        use pm_kind, only: IKC => IK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure iseq_D0_D0_PSLK5
        use pm_kind, only: LKC => LK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure iseq_D0_D0_PSLK4
        use pm_kind, only: LKC => LK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure iseq_D0_D0_PSLK3
        use pm_kind, only: LKC => LK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure iseq_D0_D0_PSLK2
        use pm_kind, only: LKC => LK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure iseq_D0_D0_PSLK1
        use pm_kind, only: LKC => LK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure iseq_D0_D0_PSCK5
        use pm_kind, only: CKC => CK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure iseq_D0_D0_PSCK4
        use pm_kind, only: CKC => CK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure iseq_D0_D0_PSCK3
        use pm_kind, only: CKC => CK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure iseq_D0_D0_PSCK2
        use pm_kind, only: CKC => CK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure iseq_D0_D0_PSCK1
        use pm_kind, only: CKC => CK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure iseq_D0_D0_PSRK5
        use pm_kind, only: RKC => RK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure iseq_D0_D0_PSRK4
        use pm_kind, only: RKC => RK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure iseq_D0_D0_PSRK3
        use pm_kind, only: RKC => RK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure iseq_D0_D0_PSRK2
        use pm_kind, only: RKC => RK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure iseq_D0_D0_PSRK1
        use pm_kind, only: RKC => RK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PS_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef iseq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define assign_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1
    module procedure assign_D0_D0_BSSK
        use pm_kind, only: SKC => SK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  SK_ENABLED

#define IK_ENABLED 1
    module procedure assign_D0_D0_BSIK
        use pm_kind, only: IKC => IK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  IK_ENABLED

#define LK_ENABLED 1
    module procedure assign_D0_D0_BSLK
        use pm_kind, only: LKC => LK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  LK_ENABLED

#define CK_ENABLED 1
    module procedure assign_D0_D0_BSCK
        use pm_kind, only: CKC => CK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  CK_ENABLED

#define RK_ENABLED 1
    module procedure assign_D0_D0_BSRK
        use pm_kind, only: RKC => RK
#include "pm_container@routines.inc.F90"
    end procedure
#undef  RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef BS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure assign_D0_D0_PSSK5
        use pm_kind, only: SKC => SK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure assign_D0_D0_PSSK4
        use pm_kind, only: SKC => SK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure assign_D0_D0_PSSK3
        use pm_kind, only: SKC => SK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure assign_D0_D0_PSSK2
        use pm_kind, only: SKC => SK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure assign_D0_D0_PSSK1
        use pm_kind, only: SKC => SK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure assign_D0_D0_PSIK5
        use pm_kind, only: IKC => IK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure assign_D0_D0_PSIK4
        use pm_kind, only: IKC => IK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure assign_D0_D0_PSIK3
        use pm_kind, only: IKC => IK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure assign_D0_D0_PSIK2
        use pm_kind, only: IKC => IK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure assign_D0_D0_PSIK1
        use pm_kind, only: IKC => IK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure assign_D0_D0_PSLK5
        use pm_kind, only: LKC => LK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure assign_D0_D0_PSLK4
        use pm_kind, only: LKC => LK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure assign_D0_D0_PSLK3
        use pm_kind, only: LKC => LK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure assign_D0_D0_PSLK2
        use pm_kind, only: LKC => LK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure assign_D0_D0_PSLK1
        use pm_kind, only: LKC => LK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure assign_D0_D0_PSCK5
        use pm_kind, only: CKC => CK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure assign_D0_D0_PSCK4
        use pm_kind, only: CKC => CK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure assign_D0_D0_PSCK3
        use pm_kind, only: CKC => CK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure assign_D0_D0_PSCK2
        use pm_kind, only: CKC => CK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure assign_D0_D0_PSCK1
        use pm_kind, only: CKC => CK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure assign_D0_D0_PSRK5
        use pm_kind, only: RKC => RK5
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure assign_D0_D0_PSRK4
        use pm_kind, only: RKC => RK4
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure assign_D0_D0_PSRK3
        use pm_kind, only: RKC => RK3
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure assign_D0_D0_PSRK2
        use pm_kind, only: RKC => RK2
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure assign_D0_D0_PSRK1
        use pm_kind, only: RKC => RK1
#include "pm_container@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PS_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef assign_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getVal_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define  D0_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!    module procedure getVal_D0_BSSK
!        use pm_kind, only: SKC => SK
!#include "pm_container@routines.inc.F90"
!    end procedure
!#undef  SK_ENABLED
!
!#define IK_ENABLED 1
!    module procedure getVal_D0_BSIK
!        use pm_kind, only: IKC => IK
!#include "pm_container@routines.inc.F90"
!    end procedure
!#undef  IK_ENABLED
!
!#define LK_ENABLED 1
!    module procedure getVal_D0_BSLK
!        use pm_kind, only: LKC => LK
!#include "pm_container@routines.inc.F90"
!    end procedure
!#undef  LK_ENABLED
!
!#define CK_ENABLED 1
!    module procedure getVal_D0_BSCK
!        use pm_kind, only: CKC => CK
!#include "pm_container@routines.inc.F90"
!    end procedure
!#undef  CK_ENABLED
!
!#define RK_ENABLED 1
!    module procedure getVal_D0_BSRK
!        use pm_kind, only: RKC => RK
!#include "pm_container@routines.inc.F90"
!    end procedure
!#undef  RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef D0_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef getVal_ENABLED
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines