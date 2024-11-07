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
!>  This file contains procedure implementations of [pm_cosmology](@ref pm_cosmology).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Wednesday 5:43 PM, December 25, 2013, Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_cosmology) routines

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define EPS 10 * epsilon(0._RKG)
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_quadPack, only: getQuadErr, GK21, pinf, weps
    use pm_val2str, only: getStr
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getSizeUnivNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Z_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getSizeUnivNormedZ_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getSizeUnivNormedZ_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getSizeUnivNormedZ_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getSizeUnivNormedZ_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getSizeUnivNormedZ_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Z_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZML_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getSizeUnivNormedZML_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getSizeUnivNormedZML_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getSizeUnivNormedZML_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getSizeUnivNormedZML_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getSizeUnivNormedZML_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZML_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZMLR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getSizeUnivNormedZMLR_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getSizeUnivNormedZMLR_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getSizeUnivNormedZMLR_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getSizeUnivNormedZMLR_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getSizeUnivNormedZMLR_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZMLR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZMLRK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getSizeUnivNormedZMLRK_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getSizeUnivNormedZMLRK_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getSizeUnivNormedZMLRK_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getSizeUnivNormedZMLRK_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getSizeUnivNormedZMLRK_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZMLRK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getSizeUnivNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisLookbackNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Z_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisLookbackNormedZ_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisLookbackNormedZ_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisLookbackNormedZ_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisLookbackNormedZ_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisLookbackNormedZ_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Z_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZML_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisLookbackNormedZML_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisLookbackNormedZML_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisLookbackNormedZML_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisLookbackNormedZML_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisLookbackNormedZML_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZML_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZMLR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisLookbackNormedZMLR_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisLookbackNormedZMLR_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisLookbackNormedZMLR_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisLookbackNormedZMLR_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisLookbackNormedZMLR_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZMLR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZMLRK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisLookbackNormedZMLRK_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisLookbackNormedZMLRK_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisLookbackNormedZMLRK_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisLookbackNormedZMLRK_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisLookbackNormedZMLRK_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZMLRK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisLookbackNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisComNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Z_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisComNormedZ_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisComNormedZ_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisComNormedZ_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisComNormedZ_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisComNormedZ_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Z_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZML_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisComNormedZML_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisComNormedZML_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisComNormedZML_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisComNormedZML_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisComNormedZML_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZML_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZMLR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisComNormedZMLR_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisComNormedZMLR_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisComNormedZMLR_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisComNormedZMLR_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisComNormedZMLR_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZMLR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZMLRK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisComNormedZMLRK_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisComNormedZMLRK_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisComNormedZMLRK_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisComNormedZMLRK_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisComNormedZMLRK_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZMLRK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisComNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisComTransNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Z_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisComTransNormedZ_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisComTransNormedZ_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisComTransNormedZ_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisComTransNormedZ_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisComTransNormedZ_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Z_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZML_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisComTransNormedZML_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisComTransNormedZML_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisComTransNormedZML_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisComTransNormedZML_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisComTransNormedZML_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZML_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZMLR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisComTransNormedZMLR_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisComTransNormedZMLR_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisComTransNormedZMLR_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisComTransNormedZMLR_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisComTransNormedZMLR_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZMLR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZMLRK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisComTransNormedZMLRK_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisComTransNormedZMLRK_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisComTransNormedZMLRK_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisComTransNormedZMLRK_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisComTransNormedZMLRK_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZMLRK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisComTransNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisAngNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Z_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisAngNormedZ_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisAngNormedZ_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisAngNormedZ_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisAngNormedZ_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisAngNormedZ_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Z_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZML_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisAngNormedZML_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisAngNormedZML_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisAngNormedZML_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisAngNormedZML_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisAngNormedZML_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZML_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZMLR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisAngNormedZMLR_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisAngNormedZMLR_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisAngNormedZMLR_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisAngNormedZMLR_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisAngNormedZMLR_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZMLR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZMLRK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisAngNormedZMLRK_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisAngNormedZMLRK_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisAngNormedZMLRK_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisAngNormedZMLRK_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisAngNormedZMLRK_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZMLRK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisAngNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisLumNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Z_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisLumNormedZ_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisLumNormedZ_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisLumNormedZ_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisLumNormedZ_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisLumNormedZ_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Z_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZML_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisLumNormedZML_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisLumNormedZML_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisLumNormedZML_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisLumNormedZML_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisLumNormedZML_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZML_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZMLR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisLumNormedZMLR_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisLumNormedZMLR_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisLumNormedZMLR_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisLumNormedZMLR_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisLumNormedZMLR_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZMLR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZMLRK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisLumNormedZMLRK_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisLumNormedZMLRK_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisLumNormedZMLRK_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisLumNormedZMLRK_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisLumNormedZMLRK_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZMLRK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisLumNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisComTransNormedWU10_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Z_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisComTransNormedWU10Z_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisComTransNormedWU10Z_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisComTransNormedWU10Z_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisComTransNormedWU10Z_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisComTransNormedWU10Z_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Z_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZML_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisComTransNormedWU10ZML_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisComTransNormedWU10ZML_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisComTransNormedWU10ZML_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisComTransNormedWU10ZML_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisComTransNormedWU10ZML_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZML_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisComTransNormedWU10_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getHubbleParamNormedSq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Z_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getHubbleParamNormedSqZ_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getHubbleParamNormedSqZ_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getHubbleParamNormedSqZ_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getHubbleParamNormedSqZ_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getHubbleParamNormedSqZ_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Z_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZML_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getHubbleParamNormedSqZML_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getHubbleParamNormedSqZML_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getHubbleParamNormedSqZML_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getHubbleParamNormedSqZML_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getHubbleParamNormedSqZML_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZML_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZMLR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getHubbleParamNormedSqZMLR_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getHubbleParamNormedSqZMLR_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getHubbleParamNormedSqZMLR_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getHubbleParamNormedSqZMLR_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getHubbleParamNormedSqZMLR_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZMLR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZMLRK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getHubbleParamNormedSqZMLRK_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getHubbleParamNormedSqZMLRK_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getHubbleParamNormedSqZMLRK_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getHubbleParamNormedSqZMLRK_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getHubbleParamNormedSqZMLRK_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZMLRK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getHubbleParamNormedSq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getVolComDiffNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Z_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVolComDiffNormedZ_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVolComDiffNormedZ_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVolComDiffNormedZ_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVolComDiffNormedZ_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVolComDiffNormedZ_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Z_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZML_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVolComDiffNormedZML_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVolComDiffNormedZML_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVolComDiffNormedZML_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVolComDiffNormedZML_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVolComDiffNormedZML_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZML_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZMLR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVolComDiffNormedZMLR_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVolComDiffNormedZMLR_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVolComDiffNormedZMLR_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVolComDiffNormedZMLR_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVolComDiffNormedZMLR_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZMLR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZMLRK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVolComDiffNormedZMLRK_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVolComDiffNormedZMLRK_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVolComDiffNormedZMLRK_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVolComDiffNormedZMLRK_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVolComDiffNormedZMLRK_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZMLRK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getVolComDiffNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setVolComDiffNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVolComDiffNormed_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVolComDiffNormed_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVolComDiffNormed_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVolComDiffNormed_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVolComDiffNormed_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setVolComDiffNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getVolComNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVolComNormedD_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVolComNormedD_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVolComNormedD_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVolComNormedD_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVolComNormedD_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DOS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVolComNormedDOS_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVolComNormedDOS_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVolComNormedDOS_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVolComNormedDOS_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVolComNormedDOS_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DOS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getVolComNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION
#undef EPS

end submodule routines