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
!>  This file contains procedure implementations of [pm_clusKmeans](@ref pm_clusKmeans).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Wednesday 5:43 PM, December 25, 2013, Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_clusKmeans) routines

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_distUnif, only: setUnifRand
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCenter_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCenterEuc_D1_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCenterEuc_D1_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCenterEuc_D1_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCenterEuc_D1_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCenterEuc_D1_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCenterEuc_D2_D2_RK5
        use pm_kind, only: RKG => RK5
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCenterEuc_D2_D2_RK4
        use pm_kind, only: RKG => RK4
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCenterEuc_D2_D2_RK3
        use pm_kind, only: RKG => RK3
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCenterEuc_D2_D2_RK2
        use pm_kind, only: RKG => RK2
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCenterEuc_D2_D2_RK1
        use pm_kind, only: RKG => RK1
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCenter_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMember_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMemberEucDef_D0_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMemberEucDef_D0_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMemberEucDef_D0_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMemberEucDef_D0_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMemberEucDef_D0_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMemberEucDef_D1_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMemberEucDef_D1_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMemberEucDef_D1_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMemberEucDef_D1_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMemberEucDef_D1_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMemberEucDef_D1_D2_RK5
        use pm_kind, only: RKG => RK5
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMemberEucDef_D1_D2_RK4
        use pm_kind, only: RKG => RK4
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMemberEucDef_D1_D2_RK3
        use pm_kind, only: RKG => RK3
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMemberEucDef_D1_D2_RK2
        use pm_kind, only: RKG => RK2
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMemberEucDef_D1_D2_RK1
        use pm_kind, only: RKG => RK1
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMemberEucDef_D2_D2_RK5
        use pm_kind, only: RKG => RK5
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMemberEucDef_D2_D2_RK4
        use pm_kind, only: RKG => RK4
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMemberEucDef_D2_D2_RK3
        use pm_kind, only: RKG => RK3
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMemberEucDef_D2_D2_RK2
        use pm_kind, only: RKG => RK2
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMemberEucDef_D2_D2_RK1
        use pm_kind, only: RKG => RK1
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Cng_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMemberEucCng_D0_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMemberEucCng_D0_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMemberEucCng_D0_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMemberEucCng_D0_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMemberEucCng_D0_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMemberEucCng_D1_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMemberEucCng_D1_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMemberEucCng_D1_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMemberEucCng_D1_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMemberEucCng_D1_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMemberEucCng_D1_D2_RK5
        use pm_kind, only: RKG => RK5
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMemberEucCng_D1_D2_RK4
        use pm_kind, only: RKG => RK4
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMemberEucCng_D1_D2_RK3
        use pm_kind, only: RKG => RK3
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMemberEucCng_D1_D2_RK2
        use pm_kind, only: RKG => RK2
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMemberEucCng_D1_D2_RK1
        use pm_kind, only: RKG => RK1
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMemberEucCng_D2_D2_RK5
        use pm_kind, only: RKG => RK5
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMemberEucCng_D2_D2_RK4
        use pm_kind, only: RKG => RK4
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMemberEucCng_D2_D2_RK3
        use pm_kind, only: RKG => RK3
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMemberEucCng_D2_D2_RK2
        use pm_kind, only: RKG => RK2
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMemberEucCng_D2_D2_RK1
        use pm_kind, only: RKG => RK1
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Cng_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMember_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setKmeansPP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Default_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setKmeansPPDRNGF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setKmeansPPDRNGF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setKmeansPPDRNGF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setKmeansPPDRNGF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setKmeansPPDRNGF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setKmeansPPDRNGX_RK5
        use pm_kind, only: RKG => RK5
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setKmeansPPDRNGX_RK4
        use pm_kind, only: RKG => RK4
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setKmeansPPDRNGX_RK3
        use pm_kind, only: RKG => RK3
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setKmeansPPDRNGX_RK2
        use pm_kind, only: RKG => RK2
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setKmeansPPDRNGX_RK1
        use pm_kind, only: RKG => RK1
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Default_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Optional_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setKmeansPPORNGF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setKmeansPPORNGF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setKmeansPPORNGF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setKmeansPPORNGF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setKmeansPPORNGF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setKmeansPPORNGX_RK5
        use pm_kind, only: RKG => RK5
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setKmeansPPORNGX_RK4
        use pm_kind, only: RKG => RK4
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setKmeansPPORNGX_RK3
        use pm_kind, only: RKG => RK3
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setKmeansPPORNGX_RK2
        use pm_kind, only: RKG => RK2
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setKmeansPPORNGX_RK1
        use pm_kind, only: RKG => RK1
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Optional_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setKmeansPP_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setKmeans_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Init_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setKmeansInit_RK5
        use pm_kind, only: RKG => RK5
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setKmeansInit_RK4
        use pm_kind, only: RKG => RK4
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setKmeansInit_RK3
        use pm_kind, only: RKG => RK3
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setKmeansInit_RK2
        use pm_kind, only: RKG => RK2
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setKmeansInit_RK1
        use pm_kind, only: RKG => RK1
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Init_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setKmeansRNGF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setKmeansRNGF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setKmeansRNGF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setKmeansRNGF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setKmeansRNGF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setKmeansRNGX_RK5
        use pm_kind, only: RKG => RK5
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setKmeansRNGX_RK4
        use pm_kind, only: RKG => RK4
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setKmeansRNGX_RK3
        use pm_kind, only: RKG => RK3
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setKmeansRNGX_RK2
        use pm_kind, only: RKG => RK2
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setKmeansRNGX_RK1
        use pm_kind, only: RKG => RK1
#include "pm_clusKmeans@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setKmeans_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines