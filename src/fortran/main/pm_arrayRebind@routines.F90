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
!>  This file contains procedure implementations of [pm_arrayRebind](@ref pm_arrayRebind).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_arrayRebind) routines ! LCOV_EXCL_LINE

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

#define setRebound_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SDDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setReboundSDDD_D1_SK5
        use pm_kind, only: SKG => SK5
        character(len(Array,IK),SKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setReboundSDDD_D1_SK4
        use pm_kind, only: SKG => SK4
        character(len(Array,IK),SKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setReboundSDDD_D1_SK3
        use pm_kind, only: SKG => SK3
        character(len(Array,IK),SKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setReboundSDDD_D1_SK2
        use pm_kind, only: SKG => SK2
        character(len(Array,IK),SKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setReboundSDDD_D1_SK1
        use pm_kind, only: SKG => SK1
        character(len(Array,IK),SKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setReboundSDDD_D1_IK5
        use pm_kind, only: IKG => IK5
        integer(IKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setReboundSDDD_D1_IK4
        use pm_kind, only: IKG => IK4
        integer(IKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setReboundSDDD_D1_IK3
        use pm_kind, only: IKG => IK3
        integer(IKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setReboundSDDD_D1_IK2
        use pm_kind, only: IKG => IK2
        integer(IKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setReboundSDDD_D1_IK1
        use pm_kind, only: IKG => IK1
        integer(IKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setReboundSDDD_D1_LK5
        use pm_kind, only: LKG => LK5
        logical(LKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setReboundSDDD_D1_LK4
        use pm_kind, only: LKG => LK4
        logical(LKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setReboundSDDD_D1_LK3
        use pm_kind, only: LKG => LK3
        logical(LKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setReboundSDDD_D1_LK2
        use pm_kind, only: LKG => LK2
        logical(LKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setReboundSDDD_D1_LK1
        use pm_kind, only: LKG => LK1
        logical(LKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setReboundSDDD_D1_CK5
        use pm_kind, only: CKG => CK5
        complex(CKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setReboundSDDD_D1_CK4
        use pm_kind, only: CKG => CK4
        complex(CKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setReboundSDDD_D1_CK3
        use pm_kind, only: CKG => CK3
        complex(CKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setReboundSDDD_D1_CK2
        use pm_kind, only: CKG => CK2
        complex(CKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setReboundSDDD_D1_CK1
        use pm_kind, only: CKG => CK1
        complex(CKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setReboundSDDD_D1_RK5
        use pm_kind, only: RKG => RK5
        real(RKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setReboundSDDD_D1_RK4
        use pm_kind, only: RKG => RK4
        real(RKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setReboundSDDD_D1_RK3
        use pm_kind, only: RKG => RK3
        real(RKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setReboundSDDD_D1_RK2
        use pm_kind, only: RKG => RK2
        real(RKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setReboundSDDD_D1_RK1
        use pm_kind, only: RKG => RK1
        real(RKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setReboundSDDD_D2_SK5
        use pm_kind, only: SKG => SK5
        character(len(Array,IK),SKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setReboundSDDD_D2_SK4
        use pm_kind, only: SKG => SK4
        character(len(Array,IK),SKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setReboundSDDD_D2_SK3
        use pm_kind, only: SKG => SK3
        character(len(Array,IK),SKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setReboundSDDD_D2_SK2
        use pm_kind, only: SKG => SK2
        character(len(Array,IK),SKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setReboundSDDD_D2_SK1
        use pm_kind, only: SKG => SK1
        character(len(Array,IK),SKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setReboundSDDD_D2_IK5
        use pm_kind, only: IKG => IK5
        integer(IKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setReboundSDDD_D2_IK4
        use pm_kind, only: IKG => IK4
        integer(IKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setReboundSDDD_D2_IK3
        use pm_kind, only: IKG => IK3
        integer(IKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setReboundSDDD_D2_IK2
        use pm_kind, only: IKG => IK2
        integer(IKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setReboundSDDD_D2_IK1
        use pm_kind, only: IKG => IK1
        integer(IKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setReboundSDDD_D2_LK5
        use pm_kind, only: LKG => LK5
        logical(LKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setReboundSDDD_D2_LK4
        use pm_kind, only: LKG => LK4
        logical(LKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setReboundSDDD_D2_LK3
        use pm_kind, only: LKG => LK3
        logical(LKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setReboundSDDD_D2_LK2
        use pm_kind, only: LKG => LK2
        logical(LKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setReboundSDDD_D2_LK1
        use pm_kind, only: LKG => LK1
        logical(LKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setReboundSDDD_D2_CK5
        use pm_kind, only: CKG => CK5
        complex(CKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setReboundSDDD_D2_CK4
        use pm_kind, only: CKG => CK4
        complex(CKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setReboundSDDD_D2_CK3
        use pm_kind, only: CKG => CK3
        complex(CKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setReboundSDDD_D2_CK2
        use pm_kind, only: CKG => CK2
        complex(CKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setReboundSDDD_D2_CK1
        use pm_kind, only: CKG => CK1
        complex(CKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setReboundSDDD_D2_RK5
        use pm_kind, only: RKG => RK5
        real(RKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setReboundSDDD_D2_RK4
        use pm_kind, only: RKG => RK4
        real(RKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setReboundSDDD_D2_RK3
        use pm_kind, only: RKG => RK3
        real(RKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setReboundSDDD_D2_RK2
        use pm_kind, only: RKG => RK2
        real(RKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setReboundSDDD_D2_RK1
        use pm_kind, only: RKG => RK1
        real(RKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setReboundSDDD_D3_SK5
        use pm_kind, only: SKG => SK5
        character(len(Array,IK),SKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setReboundSDDD_D3_SK4
        use pm_kind, only: SKG => SK4
        character(len(Array,IK),SKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setReboundSDDD_D3_SK3
        use pm_kind, only: SKG => SK3
        character(len(Array,IK),SKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setReboundSDDD_D3_SK2
        use pm_kind, only: SKG => SK2
        character(len(Array,IK),SKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setReboundSDDD_D3_SK1
        use pm_kind, only: SKG => SK1
        character(len(Array,IK),SKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setReboundSDDD_D3_IK5
        use pm_kind, only: IKG => IK5
        integer(IKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setReboundSDDD_D3_IK4
        use pm_kind, only: IKG => IK4
        integer(IKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setReboundSDDD_D3_IK3
        use pm_kind, only: IKG => IK3
        integer(IKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setReboundSDDD_D3_IK2
        use pm_kind, only: IKG => IK2
        integer(IKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setReboundSDDD_D3_IK1
        use pm_kind, only: IKG => IK1
        integer(IKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setReboundSDDD_D3_LK5
        use pm_kind, only: LKG => LK5
        logical(LKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setReboundSDDD_D3_LK4
        use pm_kind, only: LKG => LK4
        logical(LKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setReboundSDDD_D3_LK3
        use pm_kind, only: LKG => LK3
        logical(LKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setReboundSDDD_D3_LK2
        use pm_kind, only: LKG => LK2
        logical(LKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setReboundSDDD_D3_LK1
        use pm_kind, only: LKG => LK1
        logical(LKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setReboundSDDD_D3_CK5
        use pm_kind, only: CKG => CK5
        complex(CKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setReboundSDDD_D3_CK4
        use pm_kind, only: CKG => CK4
        complex(CKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setReboundSDDD_D3_CK3
        use pm_kind, only: CKG => CK3
        complex(CKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setReboundSDDD_D3_CK2
        use pm_kind, only: CKG => CK2
        complex(CKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setReboundSDDD_D3_CK1
        use pm_kind, only: CKG => CK1
        complex(CKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setReboundSDDD_D3_RK5
        use pm_kind, only: RKG => RK5
        real(RKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setReboundSDDD_D3_RK4
        use pm_kind, only: RKG => RK4
        real(RKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setReboundSDDD_D3_RK3
        use pm_kind, only: RKG => RK3
        real(RKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setReboundSDDD_D3_RK2
        use pm_kind, only: RKG => RK2
        real(RKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setReboundSDDD_D3_RK1
        use pm_kind, only: RKG => RK1
        real(RKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D3_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SDDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setReboundSLDD_D1_SK5
        use pm_kind, only: SKG => SK5
        character(len(Array,IK),SKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setReboundSLDD_D1_SK4
        use pm_kind, only: SKG => SK4
        character(len(Array,IK),SKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setReboundSLDD_D1_SK3
        use pm_kind, only: SKG => SK3
        character(len(Array,IK),SKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setReboundSLDD_D1_SK2
        use pm_kind, only: SKG => SK2
        character(len(Array,IK),SKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setReboundSLDD_D1_SK1
        use pm_kind, only: SKG => SK1
        character(len(Array,IK),SKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setReboundSLDD_D1_IK5
        use pm_kind, only: IKG => IK5
        integer(IKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setReboundSLDD_D1_IK4
        use pm_kind, only: IKG => IK4
        integer(IKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setReboundSLDD_D1_IK3
        use pm_kind, only: IKG => IK3
        integer(IKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setReboundSLDD_D1_IK2
        use pm_kind, only: IKG => IK2
        integer(IKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setReboundSLDD_D1_IK1
        use pm_kind, only: IKG => IK1
        integer(IKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setReboundSLDD_D1_LK5
        use pm_kind, only: LKG => LK5
        logical(LKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setReboundSLDD_D1_LK4
        use pm_kind, only: LKG => LK4
        logical(LKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setReboundSLDD_D1_LK3
        use pm_kind, only: LKG => LK3
        logical(LKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setReboundSLDD_D1_LK2
        use pm_kind, only: LKG => LK2
        logical(LKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setReboundSLDD_D1_LK1
        use pm_kind, only: LKG => LK1
        logical(LKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setReboundSLDD_D1_CK5
        use pm_kind, only: CKG => CK5
        complex(CKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setReboundSLDD_D1_CK4
        use pm_kind, only: CKG => CK4
        complex(CKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setReboundSLDD_D1_CK3
        use pm_kind, only: CKG => CK3
        complex(CKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setReboundSLDD_D1_CK2
        use pm_kind, only: CKG => CK2
        complex(CKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setReboundSLDD_D1_CK1
        use pm_kind, only: CKG => CK1
        complex(CKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setReboundSLDD_D1_RK5
        use pm_kind, only: RKG => RK5
        real(RKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setReboundSLDD_D1_RK4
        use pm_kind, only: RKG => RK4
        real(RKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setReboundSLDD_D1_RK3
        use pm_kind, only: RKG => RK3
        real(RKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setReboundSLDD_D1_RK2
        use pm_kind, only: RKG => RK2
        real(RKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setReboundSLDD_D1_RK1
        use pm_kind, only: RKG => RK1
        real(RKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setReboundSLDD_D2_SK5
        use pm_kind, only: SKG => SK5
        character(len(Array,IK),SKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setReboundSLDD_D2_SK4
        use pm_kind, only: SKG => SK4
        character(len(Array,IK),SKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setReboundSLDD_D2_SK3
        use pm_kind, only: SKG => SK3
        character(len(Array,IK),SKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setReboundSLDD_D2_SK2
        use pm_kind, only: SKG => SK2
        character(len(Array,IK),SKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setReboundSLDD_D2_SK1
        use pm_kind, only: SKG => SK1
        character(len(Array,IK),SKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setReboundSLDD_D2_IK5
        use pm_kind, only: IKG => IK5
        integer(IKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setReboundSLDD_D2_IK4
        use pm_kind, only: IKG => IK4
        integer(IKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setReboundSLDD_D2_IK3
        use pm_kind, only: IKG => IK3
        integer(IKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setReboundSLDD_D2_IK2
        use pm_kind, only: IKG => IK2
        integer(IKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setReboundSLDD_D2_IK1
        use pm_kind, only: IKG => IK1
        integer(IKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setReboundSLDD_D2_LK5
        use pm_kind, only: LKG => LK5
        logical(LKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setReboundSLDD_D2_LK4
        use pm_kind, only: LKG => LK4
        logical(LKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setReboundSLDD_D2_LK3
        use pm_kind, only: LKG => LK3
        logical(LKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setReboundSLDD_D2_LK2
        use pm_kind, only: LKG => LK2
        logical(LKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setReboundSLDD_D2_LK1
        use pm_kind, only: LKG => LK1
        logical(LKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setReboundSLDD_D2_CK5
        use pm_kind, only: CKG => CK5
        complex(CKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setReboundSLDD_D2_CK4
        use pm_kind, only: CKG => CK4
        complex(CKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setReboundSLDD_D2_CK3
        use pm_kind, only: CKG => CK3
        complex(CKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setReboundSLDD_D2_CK2
        use pm_kind, only: CKG => CK2
        complex(CKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setReboundSLDD_D2_CK1
        use pm_kind, only: CKG => CK1
        complex(CKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setReboundSLDD_D2_RK5
        use pm_kind, only: RKG => RK5
        real(RKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setReboundSLDD_D2_RK4
        use pm_kind, only: RKG => RK4
        real(RKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setReboundSLDD_D2_RK3
        use pm_kind, only: RKG => RK3
        real(RKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setReboundSLDD_D2_RK2
        use pm_kind, only: RKG => RK2
        real(RKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setReboundSLDD_D2_RK1
        use pm_kind, only: RKG => RK1
        real(RKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setReboundSLDD_D3_SK5
        use pm_kind, only: SKG => SK5
        character(len(Array,IK),SKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setReboundSLDD_D3_SK4
        use pm_kind, only: SKG => SK4
        character(len(Array,IK),SKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setReboundSLDD_D3_SK3
        use pm_kind, only: SKG => SK3
        character(len(Array,IK),SKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setReboundSLDD_D3_SK2
        use pm_kind, only: SKG => SK2
        character(len(Array,IK),SKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setReboundSLDD_D3_SK1
        use pm_kind, only: SKG => SK1
        character(len(Array,IK),SKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setReboundSLDD_D3_IK5
        use pm_kind, only: IKG => IK5
        integer(IKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setReboundSLDD_D3_IK4
        use pm_kind, only: IKG => IK4
        integer(IKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setReboundSLDD_D3_IK3
        use pm_kind, only: IKG => IK3
        integer(IKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setReboundSLDD_D3_IK2
        use pm_kind, only: IKG => IK2
        integer(IKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setReboundSLDD_D3_IK1
        use pm_kind, only: IKG => IK1
        integer(IKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setReboundSLDD_D3_LK5
        use pm_kind, only: LKG => LK5
        logical(LKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setReboundSLDD_D3_LK4
        use pm_kind, only: LKG => LK4
        logical(LKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setReboundSLDD_D3_LK3
        use pm_kind, only: LKG => LK3
        logical(LKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setReboundSLDD_D3_LK2
        use pm_kind, only: LKG => LK2
        logical(LKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setReboundSLDD_D3_LK1
        use pm_kind, only: LKG => LK1
        logical(LKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setReboundSLDD_D3_CK5
        use pm_kind, only: CKG => CK5
        complex(CKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setReboundSLDD_D3_CK4
        use pm_kind, only: CKG => CK4
        complex(CKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setReboundSLDD_D3_CK3
        use pm_kind, only: CKG => CK3
        complex(CKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setReboundSLDD_D3_CK2
        use pm_kind, only: CKG => CK2
        complex(CKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setReboundSLDD_D3_CK1
        use pm_kind, only: CKG => CK1
        complex(CKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setReboundSLDD_D3_RK5
        use pm_kind, only: RKG => RK5
        real(RKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setReboundSLDD_D3_RK4
        use pm_kind, only: RKG => RK4
        real(RKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setReboundSLDD_D3_RK3
        use pm_kind, only: RKG => RK3
        real(RKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setReboundSLDD_D3_RK2
        use pm_kind, only: RKG => RK2
        real(RKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setReboundSLDD_D3_RK1
        use pm_kind, only: RKG => RK1
        real(RKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D3_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SLLU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setReboundSLLU_D1_SK5
        use pm_kind, only: SKG => SK5
        character(len(Array,IK),SKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setReboundSLLU_D1_SK4
        use pm_kind, only: SKG => SK4
        character(len(Array,IK),SKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setReboundSLLU_D1_SK3
        use pm_kind, only: SKG => SK3
        character(len(Array,IK),SKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setReboundSLLU_D1_SK2
        use pm_kind, only: SKG => SK2
        character(len(Array,IK),SKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setReboundSLLU_D1_SK1
        use pm_kind, only: SKG => SK1
        character(len(Array,IK),SKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setReboundSLLU_D1_IK5
        use pm_kind, only: IKG => IK5
        integer(IKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setReboundSLLU_D1_IK4
        use pm_kind, only: IKG => IK4
        integer(IKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setReboundSLLU_D1_IK3
        use pm_kind, only: IKG => IK3
        integer(IKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setReboundSLLU_D1_IK2
        use pm_kind, only: IKG => IK2
        integer(IKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setReboundSLLU_D1_IK1
        use pm_kind, only: IKG => IK1
        integer(IKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setReboundSLLU_D1_LK5
        use pm_kind, only: LKG => LK5
        logical(LKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setReboundSLLU_D1_LK4
        use pm_kind, only: LKG => LK4
        logical(LKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setReboundSLLU_D1_LK3
        use pm_kind, only: LKG => LK3
        logical(LKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setReboundSLLU_D1_LK2
        use pm_kind, only: LKG => LK2
        logical(LKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setReboundSLLU_D1_LK1
        use pm_kind, only: LKG => LK1
        logical(LKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setReboundSLLU_D1_CK5
        use pm_kind, only: CKG => CK5
        complex(CKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setReboundSLLU_D1_CK4
        use pm_kind, only: CKG => CK4
        complex(CKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setReboundSLLU_D1_CK3
        use pm_kind, only: CKG => CK3
        complex(CKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setReboundSLLU_D1_CK2
        use pm_kind, only: CKG => CK2
        complex(CKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setReboundSLLU_D1_CK1
        use pm_kind, only: CKG => CK1
        complex(CKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setReboundSLLU_D1_RK5
        use pm_kind, only: RKG => RK5
        real(RKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setReboundSLLU_D1_RK4
        use pm_kind, only: RKG => RK4
        real(RKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setReboundSLLU_D1_RK3
        use pm_kind, only: RKG => RK3
        real(RKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setReboundSLLU_D1_RK2
        use pm_kind, only: RKG => RK2
        real(RKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setReboundSLLU_D1_RK1
        use pm_kind, only: RKG => RK1
        real(RKG), allocatable :: Temp(:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setReboundSLLU_D2_SK5
        use pm_kind, only: SKG => SK5
        character(len(Array,IK),SKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setReboundSLLU_D2_SK4
        use pm_kind, only: SKG => SK4
        character(len(Array,IK),SKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setReboundSLLU_D2_SK3
        use pm_kind, only: SKG => SK3
        character(len(Array,IK),SKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setReboundSLLU_D2_SK2
        use pm_kind, only: SKG => SK2
        character(len(Array,IK),SKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setReboundSLLU_D2_SK1
        use pm_kind, only: SKG => SK1
        character(len(Array,IK),SKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setReboundSLLU_D2_IK5
        use pm_kind, only: IKG => IK5
        integer(IKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setReboundSLLU_D2_IK4
        use pm_kind, only: IKG => IK4
        integer(IKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setReboundSLLU_D2_IK3
        use pm_kind, only: IKG => IK3
        integer(IKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setReboundSLLU_D2_IK2
        use pm_kind, only: IKG => IK2
        integer(IKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setReboundSLLU_D2_IK1
        use pm_kind, only: IKG => IK1
        integer(IKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setReboundSLLU_D2_LK5
        use pm_kind, only: LKG => LK5
        logical(LKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setReboundSLLU_D2_LK4
        use pm_kind, only: LKG => LK4
        logical(LKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setReboundSLLU_D2_LK3
        use pm_kind, only: LKG => LK3
        logical(LKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setReboundSLLU_D2_LK2
        use pm_kind, only: LKG => LK2
        logical(LKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setReboundSLLU_D2_LK1
        use pm_kind, only: LKG => LK1
        logical(LKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setReboundSLLU_D2_CK5
        use pm_kind, only: CKG => CK5
        complex(CKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setReboundSLLU_D2_CK4
        use pm_kind, only: CKG => CK4
        complex(CKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setReboundSLLU_D2_CK3
        use pm_kind, only: CKG => CK3
        complex(CKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setReboundSLLU_D2_CK2
        use pm_kind, only: CKG => CK2
        complex(CKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setReboundSLLU_D2_CK1
        use pm_kind, only: CKG => CK1
        complex(CKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setReboundSLLU_D2_RK5
        use pm_kind, only: RKG => RK5
        real(RKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setReboundSLLU_D2_RK4
        use pm_kind, only: RKG => RK4
        real(RKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setReboundSLLU_D2_RK3
        use pm_kind, only: RKG => RK3
        real(RKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setReboundSLLU_D2_RK2
        use pm_kind, only: RKG => RK2
        real(RKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setReboundSLLU_D2_RK1
        use pm_kind, only: RKG => RK1
        real(RKG), allocatable :: Temp(:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setReboundSLLU_D3_SK5
        use pm_kind, only: SKG => SK5
        character(len(Array,IK),SKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setReboundSLLU_D3_SK4
        use pm_kind, only: SKG => SK4
        character(len(Array,IK),SKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setReboundSLLU_D3_SK3
        use pm_kind, only: SKG => SK3
        character(len(Array,IK),SKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setReboundSLLU_D3_SK2
        use pm_kind, only: SKG => SK2
        character(len(Array,IK),SKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setReboundSLLU_D3_SK1
        use pm_kind, only: SKG => SK1
        character(len(Array,IK),SKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setReboundSLLU_D3_IK5
        use pm_kind, only: IKG => IK5
        integer(IKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setReboundSLLU_D3_IK4
        use pm_kind, only: IKG => IK4
        integer(IKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setReboundSLLU_D3_IK3
        use pm_kind, only: IKG => IK3
        integer(IKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setReboundSLLU_D3_IK2
        use pm_kind, only: IKG => IK2
        integer(IKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setReboundSLLU_D3_IK1
        use pm_kind, only: IKG => IK1
        integer(IKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setReboundSLLU_D3_LK5
        use pm_kind, only: LKG => LK5
        logical(LKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setReboundSLLU_D3_LK4
        use pm_kind, only: LKG => LK4
        logical(LKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setReboundSLLU_D3_LK3
        use pm_kind, only: LKG => LK3
        logical(LKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setReboundSLLU_D3_LK2
        use pm_kind, only: LKG => LK2
        logical(LKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setReboundSLLU_D3_LK1
        use pm_kind, only: LKG => LK1
        logical(LKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setReboundSLLU_D3_CK5
        use pm_kind, only: CKG => CK5
        complex(CKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setReboundSLLU_D3_CK4
        use pm_kind, only: CKG => CK4
        complex(CKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setReboundSLLU_D3_CK3
        use pm_kind, only: CKG => CK3
        complex(CKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setReboundSLLU_D3_CK2
        use pm_kind, only: CKG => CK2
        complex(CKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setReboundSLLU_D3_CK1
        use pm_kind, only: CKG => CK1
        complex(CKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setReboundSLLU_D3_RK5
        use pm_kind, only: RKG => RK5
        real(RKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setReboundSLLU_D3_RK4
        use pm_kind, only: RKG => RK4
        real(RKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setReboundSLLU_D3_RK3
        use pm_kind, only: RKG => RK3
        real(RKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setReboundSLLU_D3_RK2
        use pm_kind, only: RKG => RK2
        real(RKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setReboundSLLU_D3_RK1
        use pm_kind, only: RKG => RK1
        real(RKG), allocatable :: Temp(:,:,:)
#include "pm_arrayResize@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D3_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SLLU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRebound_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines