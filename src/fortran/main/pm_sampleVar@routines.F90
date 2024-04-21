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
!>  This file contains procedure implementations of [pm_sampleVar](@ref pm_sampleVar).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 2:48 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_sampleVar) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_sampleMean, only: getMean
    use pm_complexMinMax, only: minval, maxval
    use pm_complexCompareAll, only: operator(<=)
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getVarCorrection_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Freq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVarCorrectionFreq_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVarCorrectionFreq_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVarCorrectionFreq_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVarCorrectionFreq_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVarCorrectionFreq_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Freq_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Reli_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVarCorrectionReli_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVarCorrectionReli_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVarCorrectionReli_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVarCorrectionReli_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVarCorrectionReli_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Reli_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getVarCorrection_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define getVar_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getVarALL_WNO_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getVarALL_WNO_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getVarALL_WNO_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getVarALL_WNO_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getVarALL_WNO_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVarALL_WNO_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVarALL_WNO_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVarALL_WNO_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVarALL_WNO_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVarALL_WNO_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getVarALL_WNO_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getVarALL_WNO_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getVarALL_WNO_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getVarALL_WNO_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getVarALL_WNO_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVarALL_WNO_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVarALL_WNO_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVarALL_WNO_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVarALL_WNO_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVarALL_WNO_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getVarALL_WTI_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getVarALL_WTI_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getVarALL_WTI_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getVarALL_WTI_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getVarALL_WTI_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVarALL_WTI_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVarALL_WTI_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVarALL_WTI_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVarALL_WTI_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVarALL_WTI_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getVarALL_WTI_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getVarALL_WTI_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getVarALL_WTI_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getVarALL_WTI_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getVarALL_WTI_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVarALL_WTI_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVarALL_WTI_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVarALL_WTI_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVarALL_WTI_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVarALL_WTI_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getVarALL_WTR_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getVarALL_WTR_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getVarALL_WTR_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getVarALL_WTR_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getVarALL_WTR_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVarALL_WTR_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVarALL_WTR_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVarALL_WTR_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVarALL_WTR_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVarALL_WTR_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getVarALL_WTR_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getVarALL_WTR_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getVarALL_WTR_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getVarALL_WTR_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getVarALL_WTR_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVarALL_WTR_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVarALL_WTR_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVarALL_WTR_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVarALL_WTR_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVarALL_WTR_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getVar_ENABLED
#endif
!FORTRAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define getVar_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DIM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getVarDIM_WNO_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getVarDIM_WNO_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getVarDIM_WNO_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getVarDIM_WNO_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getVarDIM_WNO_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVarDIM_WNO_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVarDIM_WNO_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVarDIM_WNO_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVarDIM_WNO_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVarDIM_WNO_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getVarDIM_WNO_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getVarDIM_WNO_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getVarDIM_WNO_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getVarDIM_WNO_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getVarDIM_WNO_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVarDIM_WNO_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVarDIM_WNO_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVarDIM_WNO_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVarDIM_WNO_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVarDIM_WNO_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getVarDIM_WTI_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getVarDIM_WTI_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getVarDIM_WTI_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getVarDIM_WTI_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getVarDIM_WTI_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVarDIM_WTI_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVarDIM_WTI_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVarDIM_WTI_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVarDIM_WTI_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVarDIM_WTI_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getVarDIM_WTI_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getVarDIM_WTI_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getVarDIM_WTI_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getVarDIM_WTI_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getVarDIM_WTI_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVarDIM_WTI_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVarDIM_WTI_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVarDIM_WTI_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVarDIM_WTI_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVarDIM_WTI_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getVarDIM_WTR_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getVarDIM_WTR_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getVarDIM_WTR_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getVarDIM_WTR_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getVarDIM_WTR_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVarDIM_WTR_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVarDIM_WTR_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVarDIM_WTR_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVarDIM_WTR_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVarDIM_WTR_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getVarDIM_WTR_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getVarDIM_WTR_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getVarDIM_WTR_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getVarDIM_WTR_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getVarDIM_WTR_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVarDIM_WTR_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVarDIM_WTR_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVarDIM_WTR_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVarDIM_WTR_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVarDIM_WTR_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DIM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getVar_ENABLED
#endif
!FORTRAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define setVar_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Avg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarAvgALL_WNO_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarAvgALL_WNO_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarAvgALL_WNO_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarAvgALL_WNO_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarAvgALL_WNO_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarAvgALL_WNO_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarAvgALL_WNO_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarAvgALL_WNO_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarAvgALL_WNO_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarAvgALL_WNO_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarAvgALL_WNO_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarAvgALL_WNO_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarAvgALL_WNO_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarAvgALL_WNO_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarAvgALL_WNO_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarAvgALL_WNO_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarAvgALL_WNO_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarAvgALL_WNO_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarAvgALL_WNO_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarAvgALL_WNO_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarAvgALL_WTI_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarAvgALL_WTI_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarAvgALL_WTI_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarAvgALL_WTI_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarAvgALL_WTI_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarAvgALL_WTI_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarAvgALL_WTI_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarAvgALL_WTI_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarAvgALL_WTI_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarAvgALL_WTI_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarAvgALL_WTI_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarAvgALL_WTI_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarAvgALL_WTI_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarAvgALL_WTI_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarAvgALL_WTI_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarAvgALL_WTI_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarAvgALL_WTI_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarAvgALL_WTI_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarAvgALL_WTI_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarAvgALL_WTI_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarAvgALL_WTR_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarAvgALL_WTR_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarAvgALL_WTR_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarAvgALL_WTR_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarAvgALL_WTR_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarAvgALL_WTR_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarAvgALL_WTR_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarAvgALL_WTR_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarAvgALL_WTR_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarAvgALL_WTR_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarAvgALL_WTR_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarAvgALL_WTR_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarAvgALL_WTR_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarAvgALL_WTR_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarAvgALL_WTR_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarAvgALL_WTR_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarAvgALL_WTR_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarAvgALL_WTR_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarAvgALL_WTR_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarAvgALL_WTR_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Avg_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Org_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarOrgALL_WNO_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarOrgALL_WNO_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarOrgALL_WNO_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarOrgALL_WNO_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarOrgALL_WNO_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarOrgALL_WNO_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarOrgALL_WNO_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarOrgALL_WNO_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarOrgALL_WNO_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarOrgALL_WNO_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarOrgALL_WNO_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarOrgALL_WNO_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarOrgALL_WNO_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarOrgALL_WNO_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarOrgALL_WNO_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarOrgALL_WNO_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarOrgALL_WNO_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarOrgALL_WNO_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarOrgALL_WNO_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarOrgALL_WNO_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarOrgALL_WTI_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarOrgALL_WTI_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarOrgALL_WTI_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarOrgALL_WTI_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarOrgALL_WTI_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarOrgALL_WTI_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarOrgALL_WTI_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarOrgALL_WTI_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarOrgALL_WTI_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarOrgALL_WTI_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarOrgALL_WTI_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarOrgALL_WTI_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarOrgALL_WTI_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarOrgALL_WTI_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarOrgALL_WTI_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarOrgALL_WTI_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarOrgALL_WTI_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarOrgALL_WTI_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarOrgALL_WTI_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarOrgALL_WTI_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarOrgALL_WTR_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarOrgALL_WTR_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarOrgALL_WTR_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarOrgALL_WTR_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarOrgALL_WTR_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarOrgALL_WTR_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarOrgALL_WTR_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarOrgALL_WTR_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarOrgALL_WTR_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarOrgALL_WTR_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarOrgALL_WTR_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarOrgALL_WTR_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarOrgALL_WTR_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarOrgALL_WTR_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarOrgALL_WTR_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarOrgALL_WTR_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarOrgALL_WTR_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarOrgALL_WTR_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarOrgALL_WTR_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarOrgALL_WTR_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Org_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setVar_ENABLED
#endif
!FORTRAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define setVar_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DIM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Avg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarAvgDIM_WNO_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarAvgDIM_WNO_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarAvgDIM_WNO_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarAvgDIM_WNO_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarAvgDIM_WNO_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarAvgDIM_WNO_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarAvgDIM_WNO_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarAvgDIM_WNO_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarAvgDIM_WNO_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarAvgDIM_WNO_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarAvgDIM_WNO_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarAvgDIM_WNO_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarAvgDIM_WNO_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarAvgDIM_WNO_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarAvgDIM_WNO_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarAvgDIM_WNO_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarAvgDIM_WNO_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarAvgDIM_WNO_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarAvgDIM_WNO_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarAvgDIM_WNO_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarAvgDIM_WTI_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarAvgDIM_WTI_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarAvgDIM_WTI_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarAvgDIM_WTI_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarAvgDIM_WTI_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarAvgDIM_WTI_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarAvgDIM_WTI_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarAvgDIM_WTI_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarAvgDIM_WTI_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarAvgDIM_WTI_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarAvgDIM_WTI_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarAvgDIM_WTI_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarAvgDIM_WTI_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarAvgDIM_WTI_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarAvgDIM_WTI_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarAvgDIM_WTI_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarAvgDIM_WTI_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarAvgDIM_WTI_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarAvgDIM_WTI_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarAvgDIM_WTI_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarAvgDIM_WTR_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarAvgDIM_WTR_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarAvgDIM_WTR_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarAvgDIM_WTR_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarAvgDIM_WTR_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarAvgDIM_WTR_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarAvgDIM_WTR_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarAvgDIM_WTR_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarAvgDIM_WTR_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarAvgDIM_WTR_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarAvgDIM_WTR_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarAvgDIM_WTR_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarAvgDIM_WTR_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarAvgDIM_WTR_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarAvgDIM_WTR_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarAvgDIM_WTR_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarAvgDIM_WTR_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarAvgDIM_WTR_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarAvgDIM_WTR_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarAvgDIM_WTR_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Avg_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Org_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarOrgDIM_WNO_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarOrgDIM_WNO_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarOrgDIM_WNO_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarOrgDIM_WNO_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarOrgDIM_WNO_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarOrgDIM_WNO_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarOrgDIM_WNO_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarOrgDIM_WNO_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarOrgDIM_WNO_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarOrgDIM_WNO_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarOrgDIM_WNO_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarOrgDIM_WNO_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarOrgDIM_WNO_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarOrgDIM_WNO_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarOrgDIM_WNO_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarOrgDIM_WNO_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarOrgDIM_WNO_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarOrgDIM_WNO_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarOrgDIM_WNO_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarOrgDIM_WNO_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarOrgDIM_WTI_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarOrgDIM_WTI_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarOrgDIM_WTI_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarOrgDIM_WTI_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarOrgDIM_WTI_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarOrgDIM_WTI_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarOrgDIM_WTI_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarOrgDIM_WTI_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarOrgDIM_WTI_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarOrgDIM_WTI_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarOrgDIM_WTI_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarOrgDIM_WTI_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarOrgDIM_WTI_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarOrgDIM_WTI_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarOrgDIM_WTI_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarOrgDIM_WTI_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarOrgDIM_WTI_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarOrgDIM_WTI_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarOrgDIM_WTI_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarOrgDIM_WTI_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarOrgDIM_WTR_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarOrgDIM_WTR_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarOrgDIM_WTR_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarOrgDIM_WTR_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarOrgDIM_WTR_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarOrgDIM_WTR_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarOrgDIM_WTR_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarOrgDIM_WTR_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarOrgDIM_WTR_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarOrgDIM_WTR_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarOrgDIM_WTR_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarOrgDIM_WTR_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarOrgDIM_WTR_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarOrgDIM_WTR_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarOrgDIM_WTR_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarOrgDIM_WTR_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarOrgDIM_WTR_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarOrgDIM_WTR_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarOrgDIM_WTR_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarOrgDIM_WTR_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Org_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DIM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setVar_ENABLED
#endif
!FORTRAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define setVarMean_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMeanALL_WNODD_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMeanALL_WNODD_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMeanALL_WNODD_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMeanALL_WNODD_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMeanALL_WNODD_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMeanALL_WNODD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMeanALL_WNODD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMeanALL_WNODD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMeanALL_WNODD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMeanALL_WNODD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMeanALL_WNODD_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMeanALL_WNODD_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMeanALL_WNODD_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMeanALL_WNODD_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMeanALL_WNODD_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMeanALL_WNODD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMeanALL_WNODD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMeanALL_WNODD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMeanALL_WNODD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMeanALL_WNODD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMeanALL_WTISD_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMeanALL_WTISD_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMeanALL_WTISD_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMeanALL_WTISD_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMeanALL_WTISD_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMeanALL_WTISD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMeanALL_WTISD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMeanALL_WTISD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMeanALL_WTISD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMeanALL_WTISD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMeanALL_WTISD_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMeanALL_WTISD_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMeanALL_WTISD_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMeanALL_WTISD_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMeanALL_WTISD_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMeanALL_WTISD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMeanALL_WTISD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMeanALL_WTISD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMeanALL_WTISD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMeanALL_WTISD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMeanALL_WTRSD_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMeanALL_WTRSD_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMeanALL_WTRSD_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMeanALL_WTRSD_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMeanALL_WTRSD_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMeanALL_WTRSD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMeanALL_WTRSD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMeanALL_WTRSD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMeanALL_WTRSD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMeanALL_WTRSD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMeanALL_WTRSD_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMeanALL_WTRSD_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMeanALL_WTRSD_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMeanALL_WTRSD_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMeanALL_WTRSD_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMeanALL_WTRSD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMeanALL_WTRSD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMeanALL_WTRSD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMeanALL_WTRSD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMeanALL_WTRSD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setVarMean_ENABLED
#endif
!FORTRAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define setVarMean_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DIM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMeanDIM_WNODD_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMeanDIM_WNODD_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMeanDIM_WNODD_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMeanDIM_WNODD_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMeanDIM_WNODD_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMeanDIM_WNODD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMeanDIM_WNODD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMeanDIM_WNODD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMeanDIM_WNODD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMeanDIM_WNODD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMeanDIM_WNODD_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMeanDIM_WNODD_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMeanDIM_WNODD_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMeanDIM_WNODD_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMeanDIM_WNODD_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMeanDIM_WNODD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMeanDIM_WNODD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMeanDIM_WNODD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMeanDIM_WNODD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMeanDIM_WNODD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMeanDIM_WTISD_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMeanDIM_WTISD_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMeanDIM_WTISD_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMeanDIM_WTISD_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMeanDIM_WTISD_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMeanDIM_WTISD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMeanDIM_WTISD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMeanDIM_WTISD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMeanDIM_WTISD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMeanDIM_WTISD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMeanDIM_WTISD_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMeanDIM_WTISD_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMeanDIM_WTISD_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMeanDIM_WTISD_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMeanDIM_WTISD_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMeanDIM_WTISD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMeanDIM_WTISD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMeanDIM_WTISD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMeanDIM_WTISD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMeanDIM_WTISD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMeanDIM_WTRSD_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMeanDIM_WTRSD_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMeanDIM_WTRSD_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMeanDIM_WTRSD_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMeanDIM_WTRSD_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMeanDIM_WTRSD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMeanDIM_WTRSD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMeanDIM_WTRSD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMeanDIM_WTRSD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMeanDIM_WTRSD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMeanDIM_WTRSD_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMeanDIM_WTRSD_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMeanDIM_WTRSD_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMeanDIM_WTRSD_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMeanDIM_WTRSD_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMeanDIM_WTRSD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMeanDIM_WTRSD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMeanDIM_WTRSD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMeanDIM_WTRSD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMeanDIM_WTRSD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DIM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setVarMean_ENABLED
#endif
!FORTRAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define getVarMerged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getVarMergedNew_D0_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getVarMergedNew_D0_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getVarMergedNew_D0_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getVarMergedNew_D0_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getVarMergedNew_D0_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVarMergedNew_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVarMergedNew_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVarMergedNew_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVarMergedNew_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVarMergedNew_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getVarMergedNew_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getVarMergedNew_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getVarMergedNew_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getVarMergedNew_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getVarMergedNew_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVarMergedNew_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVarMergedNew_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVarMergedNew_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVarMergedNew_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVarMergedNew_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getVarMerged_ENABLED
#endif
!FORTRAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define setVarMerged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMergedNew_D0_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMergedNew_D0_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMergedNew_D0_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMergedNew_D0_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMergedNew_D0_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMergedNew_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMergedNew_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMergedNew_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMergedNew_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMergedNew_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMergedNew_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMergedNew_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMergedNew_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMergedNew_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMergedNew_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMergedNew_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMergedNew_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMergedNew_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMergedNew_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMergedNew_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Old_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMergedOld_D0_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMergedOld_D0_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMergedOld_D0_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMergedOld_D0_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMergedOld_D0_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMergedOld_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMergedOld_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMergedOld_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMergedOld_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMergedOld_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMergedOld_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMergedOld_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMergedOld_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMergedOld_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMergedOld_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMergedOld_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMergedOld_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMergedOld_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMergedOld_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMergedOld_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Old_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setVarMerged_ENABLED
#endif
!FORTRAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define setVarMeanMerged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMeanMergedNew_D0_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMeanMergedNew_D0_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMeanMergedNew_D0_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMeanMergedNew_D0_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMeanMergedNew_D0_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMeanMergedNew_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMeanMergedNew_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMeanMergedNew_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMeanMergedNew_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMeanMergedNew_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMeanMergedNew_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMeanMergedNew_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMeanMergedNew_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMeanMergedNew_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMeanMergedNew_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMeanMergedNew_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMeanMergedNew_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMeanMergedNew_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMeanMergedNew_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMeanMergedNew_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Old_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMeanMergedOld_D0_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMeanMergedOld_D0_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMeanMergedOld_D0_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMeanMergedOld_D0_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMeanMergedOld_D0_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMeanMergedOld_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMeanMergedOld_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMeanMergedOld_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMeanMergedOld_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMeanMergedOld_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setVarMeanMergedOld_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setVarMeanMergedOld_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setVarMeanMergedOld_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setVarMeanMergedOld_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setVarMeanMergedOld_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVarMeanMergedOld_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVarMeanMergedOld_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVarMeanMergedOld_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVarMeanMergedOld_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVarMeanMergedOld_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Old_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setVarMeanMerged_ENABLED
#endif
!FORTRAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines