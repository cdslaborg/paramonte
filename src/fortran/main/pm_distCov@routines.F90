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
!>  This file contains procedure implementations of [pm_distCov](@ref pm_distCov).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distCov) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_array, only: nothing
    use pm_matrixChol, only: setMatChol
    use pm_distBeta, only: setBetaRand
    use pm_distUnif, only: setUnifRand
    use pm_distNorm, only: setNormRand
    use pm_matrixCopy, only: setMatCopy, upp 
    use pm_matrixCopy, only: rdpack, transHerm 
    use pm_sampleCor, only: setCor, uppDia

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCovRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GRAM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define S0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getCovRandGRNGDS0_CK5
        use pm_kind, only: TKC => CK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCovRandGRNGDS0_CK4
        use pm_kind, only: TKC => CK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCovRandGRNGDS0_CK3
        use pm_kind, only: TKC => CK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCovRandGRNGDS0_CK2
        use pm_kind, only: TKC => CK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCovRandGRNGDS0_CK1
        use pm_kind, only: TKC => CK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCovRandGRNGDS0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCovRandGRNGDS0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCovRandGRNGDS0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCovRandGRNGDS0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCovRandGRNGDS0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef S0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define S1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getCovRandGRNGDS1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCovRandGRNGDS1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCovRandGRNGDS1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCovRandGRNGDS1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCovRandGRNGDS1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCovRandGRNGDS1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCovRandGRNGDS1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCovRandGRNGDS1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCovRandGRNGDS1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCovRandGRNGDS1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef S1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GRAM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCovRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCovRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GRAM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovRandGRNGFSD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovRandGRNGFSD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovRandGRNGFSD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovRandGRNGFSD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovRandGRNGFSD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandGRNGFSD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandGRNGFSD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandGRNGFSD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandGRNGFSD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandGRNGFSD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define S0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovRandGRNGFS0_CK5
        use pm_kind, only: TKC => CK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovRandGRNGFS0_CK4
        use pm_kind, only: TKC => CK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovRandGRNGFS0_CK3
        use pm_kind, only: TKC => CK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovRandGRNGFS0_CK2
        use pm_kind, only: TKC => CK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovRandGRNGFS0_CK1
        use pm_kind, only: TKC => CK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandGRNGFS0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandGRNGFS0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandGRNGFS0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandGRNGFS0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandGRNGFS0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef S0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define S1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovRandGRNGFS1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovRandGRNGFS1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovRandGRNGFS1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovRandGRNGFS1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovRandGRNGFS1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandGRNGFS1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandGRNGFS1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandGRNGFS1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandGRNGFS1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandGRNGFS1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef S1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovRandGRNGXSD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovRandGRNGXSD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovRandGRNGXSD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovRandGRNGXSD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovRandGRNGXSD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandGRNGXSD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandGRNGXSD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandGRNGXSD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandGRNGXSD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandGRNGXSD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define S0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovRandGRNGXS0_CK5
        use pm_kind, only: TKC => CK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovRandGRNGXS0_CK4
        use pm_kind, only: TKC => CK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovRandGRNGXS0_CK3
        use pm_kind, only: TKC => CK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovRandGRNGXS0_CK2
        use pm_kind, only: TKC => CK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovRandGRNGXS0_CK1
        use pm_kind, only: TKC => CK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandGRNGXS0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandGRNGXS0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandGRNGXS0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandGRNGXS0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandGRNGXS0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef S0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define S1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovRandGRNGXS1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovRandGRNGXS1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovRandGRNGXS1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovRandGRNGXS1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovRandGRNGXS1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandGRNGXS1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandGRNGXS1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandGRNGXS1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandGRNGXS1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandGRNGXS1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef S1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GRAM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DVINE_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandDRNGFSD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandDRNGFSD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandDRNGFSD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandDRNGFSD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandDRNGFSD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define S0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandDRNGFS0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandDRNGFS0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandDRNGFS0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandDRNGFS0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandDRNGFS0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef S0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define S1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandDRNGFS1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandDRNGFS1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandDRNGFS1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandDRNGFS1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandDRNGFS1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef S1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandDRNGXSD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandDRNGXSD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandDRNGXSD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandDRNGXSD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandDRNGXSD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define S0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandDRNGXS0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandDRNGXS0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandDRNGXS0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandDRNGXS0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandDRNGXS0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef S0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define S1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandDRNGXS1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandDRNGXS1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandDRNGXS1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandDRNGXS1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandDRNGXS1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef S1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DVINE_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONION_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandORNGFSD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandORNGFSD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandORNGFSD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandORNGFSD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandORNGFSD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define S0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandORNGFS0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandORNGFS0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandORNGFS0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandORNGFS0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandORNGFS0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef S0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define S1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandORNGFS1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandORNGFS1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandORNGFS1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandORNGFS1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandORNGFS1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef S1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandORNGXSD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandORNGXSD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandORNGXSD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandORNGXSD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandORNGXSD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define S0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandORNGXS0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandORNGXS0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandORNGXS0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandORNGXS0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandORNGXS0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef S0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define S1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovRandORNGXS1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovRandORNGXS1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovRandORNGXS1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovRandORNGXS1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovRandORNGXS1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef S1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONION_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCovRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines