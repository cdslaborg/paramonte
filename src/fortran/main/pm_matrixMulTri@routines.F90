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
!>  This file contains procedure implementations of [pm_matrixMulTri](@ref pm_matrixMulTri).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_matrixMulTri) routines ! LCOV_EXCL_LINE

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
    use pm_blas, only: blasTRMV, blasTRSV, blasTRMM, blasTRSM

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatMulTri_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define trmv_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ASS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CGMB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLDA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_ASS_CLDA_ONOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_ASS_CLDA_ONOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_ASS_CLDA_ONOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_ASS_CLDA_ONOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_ASS_CLDA_ONOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_ASS_CLDA_ONOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_ASS_CLDA_ONOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_ASS_CLDA_ONOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_ASS_CLDA_ONOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_ASS_CLDA_ONOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_ASS_CLDA_OTSA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_ASS_CLDA_OTSA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_ASS_CLDA_OTSA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_ASS_CLDA_OTSA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_ASS_CLDA_OTSA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_ASS_CLDA_OTSA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_ASS_CLDA_OTSA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_ASS_CLDA_OTSA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_ASS_CLDA_OTSA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_ASS_CLDA_OTSA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_ASS_CLDA_OTHA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_ASS_CLDA_OTHA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_ASS_CLDA_OTHA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_ASS_CLDA_OTHA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_ASS_CLDA_OTHA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_ASS_CLDA_OTHA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_ASS_CLDA_OTHA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_ASS_CLDA_OTHA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_ASS_CLDA_OTHA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_ASS_CLDA_OTHA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLDA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUDA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_ASS_CUDA_ONOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_ASS_CUDA_ONOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_ASS_CUDA_ONOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_ASS_CUDA_ONOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_ASS_CUDA_ONOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_ASS_CUDA_ONOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_ASS_CUDA_ONOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_ASS_CUDA_ONOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_ASS_CUDA_ONOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_ASS_CUDA_ONOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_ASS_CUDA_OTSA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_ASS_CUDA_OTSA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_ASS_CUDA_OTSA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_ASS_CUDA_OTSA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_ASS_CUDA_OTSA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_ASS_CUDA_OTSA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_ASS_CUDA_OTSA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_ASS_CUDA_OTSA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_ASS_CUDA_OTSA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_ASS_CUDA_OTSA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_ASS_CUDA_OTHA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_ASS_CUDA_OTHA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_ASS_CUDA_OTHA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_ASS_CUDA_OTHA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_ASS_CUDA_OTHA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_ASS_CUDA_OTHA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_ASS_CUDA_OTHA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_ASS_CUDA_OTHA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_ASS_CUDA_OTHA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_ASS_CUDA_OTHA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUDA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_ASS_CLUA_ONOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_ASS_CLUA_ONOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_ASS_CLUA_ONOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_ASS_CLUA_ONOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_ASS_CLUA_ONOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_ASS_CLUA_ONOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_ASS_CLUA_ONOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_ASS_CLUA_ONOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_ASS_CLUA_ONOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_ASS_CLUA_ONOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_ASS_CLUA_OTSA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_ASS_CLUA_OTSA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_ASS_CLUA_OTSA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_ASS_CLUA_OTSA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_ASS_CLUA_OTSA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_ASS_CLUA_OTSA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_ASS_CLUA_OTSA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_ASS_CLUA_OTSA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_ASS_CLUA_OTSA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_ASS_CLUA_OTSA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_ASS_CLUA_OTHA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_ASS_CLUA_OTHA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_ASS_CLUA_OTHA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_ASS_CLUA_OTHA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_ASS_CLUA_OTHA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_ASS_CLUA_OTHA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_ASS_CLUA_OTHA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_ASS_CLUA_OTHA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_ASS_CLUA_OTHA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_ASS_CLUA_OTHA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_ASS_CUUA_ONOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_ASS_CUUA_ONOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_ASS_CUUA_ONOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_ASS_CUUA_ONOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_ASS_CUUA_ONOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_ASS_CUUA_ONOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_ASS_CUUA_ONOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_ASS_CUUA_ONOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_ASS_CUUA_ONOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_ASS_CUUA_ONOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_ASS_CUUA_OTSA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_ASS_CUUA_OTSA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_ASS_CUUA_OTSA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_ASS_CUUA_OTSA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_ASS_CUUA_OTSA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_ASS_CUUA_OTSA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_ASS_CUUA_OTSA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_ASS_CUUA_OTSA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_ASS_CUUA_OTSA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_ASS_CUUA_OTSA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_ASS_CUUA_OTHA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_ASS_CUUA_OTHA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_ASS_CUUA_OTHA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_ASS_CUUA_OTHA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_ASS_CUUA_OTHA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_ASS_CUUA_OTHA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_ASS_CUUA_OTHA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_ASS_CUUA_OTHA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_ASS_CUUA_OTHA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_ASS_CUUA_OTHA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CGMB_ENABLED

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

#define CGMB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLDA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_EXP_CLDA_ONOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_EXP_CLDA_ONOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_EXP_CLDA_ONOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_EXP_CLDA_ONOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_EXP_CLDA_ONOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_EXP_CLDA_ONOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_EXP_CLDA_ONOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_EXP_CLDA_ONOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_EXP_CLDA_ONOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_EXP_CLDA_ONOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_EXP_CLDA_OTSA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_EXP_CLDA_OTSA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_EXP_CLDA_OTSA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_EXP_CLDA_OTSA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_EXP_CLDA_OTSA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_EXP_CLDA_OTSA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_EXP_CLDA_OTSA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_EXP_CLDA_OTSA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_EXP_CLDA_OTSA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_EXP_CLDA_OTSA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_EXP_CLDA_OTHA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_EXP_CLDA_OTHA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_EXP_CLDA_OTHA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_EXP_CLDA_OTHA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_EXP_CLDA_OTHA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_EXP_CLDA_OTHA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_EXP_CLDA_OTHA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_EXP_CLDA_OTHA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_EXP_CLDA_OTHA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_EXP_CLDA_OTHA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLDA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUDA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_EXP_CUDA_ONOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_EXP_CUDA_ONOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_EXP_CUDA_ONOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_EXP_CUDA_ONOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_EXP_CUDA_ONOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_EXP_CUDA_ONOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_EXP_CUDA_ONOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_EXP_CUDA_ONOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_EXP_CUDA_ONOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_EXP_CUDA_ONOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_EXP_CUDA_OTSA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_EXP_CUDA_OTSA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_EXP_CUDA_OTSA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_EXP_CUDA_OTSA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_EXP_CUDA_OTSA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_EXP_CUDA_OTSA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_EXP_CUDA_OTSA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_EXP_CUDA_OTSA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_EXP_CUDA_OTSA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_EXP_CUDA_OTSA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_EXP_CUDA_OTHA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_EXP_CUDA_OTHA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_EXP_CUDA_OTHA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_EXP_CUDA_OTHA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_EXP_CUDA_OTHA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_EXP_CUDA_OTHA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_EXP_CUDA_OTHA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_EXP_CUDA_OTHA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_EXP_CUDA_OTHA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_EXP_CUDA_OTHA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUDA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_EXP_CLUA_ONOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_EXP_CLUA_ONOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_EXP_CLUA_ONOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_EXP_CLUA_ONOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_EXP_CLUA_ONOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_EXP_CLUA_ONOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_EXP_CLUA_ONOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_EXP_CLUA_ONOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_EXP_CLUA_ONOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_EXP_CLUA_ONOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_EXP_CLUA_OTSA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_EXP_CLUA_OTSA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_EXP_CLUA_OTSA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_EXP_CLUA_OTSA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_EXP_CLUA_OTSA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_EXP_CLUA_OTSA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_EXP_CLUA_OTSA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_EXP_CLUA_OTSA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_EXP_CLUA_OTSA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_EXP_CLUA_OTSA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_EXP_CLUA_OTHA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_EXP_CLUA_OTHA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_EXP_CLUA_OTHA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_EXP_CLUA_OTHA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_EXP_CLUA_OTHA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_EXP_CLUA_OTHA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_EXP_CLUA_OTHA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_EXP_CLUA_OTHA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_EXP_CLUA_OTHA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_EXP_CLUA_OTHA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_EXP_CUUA_ONOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_EXP_CUUA_ONOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_EXP_CUUA_ONOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_EXP_CUUA_ONOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_EXP_CUUA_ONOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_EXP_CUUA_ONOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_EXP_CUUA_ONOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_EXP_CUUA_ONOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_EXP_CUUA_ONOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_EXP_CUUA_ONOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_EXP_CUUA_OTSA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_EXP_CUUA_OTSA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_EXP_CUUA_OTSA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_EXP_CUUA_OTSA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_EXP_CUUA_OTSA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_EXP_CUUA_OTSA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_EXP_CUUA_OTSA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_EXP_CUUA_OTSA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_EXP_CUUA_OTSA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_EXP_CUUA_OTSA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmv_EXP_CUUA_OTHA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmv_EXP_CUUA_OTHA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmv_EXP_CUUA_OTHA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmv_EXP_CUUA_OTHA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmv_EXP_CUUA_OTHA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmv_EXP_CUUA_OTHA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmv_EXP_CUUA_OTHA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmv_EXP_CUUA_OTHA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmv_EXP_CUUA_OTHA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmv_EXP_CUUA_OTHA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CGMB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef EXP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef trmv_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatMulTri_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatMulTri_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define trsv_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ASS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CGMB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLDA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_ASS_CLDA_INVA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_ASS_CLDA_INVA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_ASS_CLDA_INVA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_ASS_CLDA_INVA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_ASS_CLDA_INVA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_ASS_CLDA_INVA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_ASS_CLDA_INVA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_ASS_CLDA_INVA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_ASS_CLDA_INVA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_ASS_CLDA_INVA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_ASS_CLDA_OTOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_ASS_CLDA_OTOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_ASS_CLDA_OTOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_ASS_CLDA_OTOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_ASS_CLDA_OTOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_ASS_CLDA_OTOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_ASS_CLDA_OTOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_ASS_CLDA_OTOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_ASS_CLDA_OTOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_ASS_CLDA_OTOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_ASS_CLDA_OTUA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_ASS_CLDA_OTUA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_ASS_CLDA_OTUA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_ASS_CLDA_OTUA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_ASS_CLDA_OTUA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_ASS_CLDA_OTUA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_ASS_CLDA_OTUA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_ASS_CLDA_OTUA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_ASS_CLDA_OTUA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_ASS_CLDA_OTUA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLDA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUDA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_ASS_CUDA_INVA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_ASS_CUDA_INVA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_ASS_CUDA_INVA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_ASS_CUDA_INVA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_ASS_CUDA_INVA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_ASS_CUDA_INVA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_ASS_CUDA_INVA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_ASS_CUDA_INVA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_ASS_CUDA_INVA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_ASS_CUDA_INVA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_ASS_CUDA_OTOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_ASS_CUDA_OTOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_ASS_CUDA_OTOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_ASS_CUDA_OTOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_ASS_CUDA_OTOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_ASS_CUDA_OTOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_ASS_CUDA_OTOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_ASS_CUDA_OTOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_ASS_CUDA_OTOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_ASS_CUDA_OTOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_ASS_CUDA_OTUA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_ASS_CUDA_OTUA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_ASS_CUDA_OTUA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_ASS_CUDA_OTUA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_ASS_CUDA_OTUA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_ASS_CUDA_OTUA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_ASS_CUDA_OTUA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_ASS_CUDA_OTUA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_ASS_CUDA_OTUA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_ASS_CUDA_OTUA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUDA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_ASS_CLUA_INVA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_ASS_CLUA_INVA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_ASS_CLUA_INVA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_ASS_CLUA_INVA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_ASS_CLUA_INVA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_ASS_CLUA_INVA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_ASS_CLUA_INVA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_ASS_CLUA_INVA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_ASS_CLUA_INVA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_ASS_CLUA_INVA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_ASS_CLUA_OTOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_ASS_CLUA_OTOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_ASS_CLUA_OTOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_ASS_CLUA_OTOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_ASS_CLUA_OTOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_ASS_CLUA_OTOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_ASS_CLUA_OTOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_ASS_CLUA_OTOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_ASS_CLUA_OTOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_ASS_CLUA_OTOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_ASS_CLUA_OTUA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_ASS_CLUA_OTUA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_ASS_CLUA_OTUA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_ASS_CLUA_OTUA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_ASS_CLUA_OTUA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_ASS_CLUA_OTUA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_ASS_CLUA_OTUA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_ASS_CLUA_OTUA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_ASS_CLUA_OTUA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_ASS_CLUA_OTUA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_ASS_CUUA_INVA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_ASS_CUUA_INVA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_ASS_CUUA_INVA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_ASS_CUUA_INVA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_ASS_CUUA_INVA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_ASS_CUUA_INVA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_ASS_CUUA_INVA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_ASS_CUUA_INVA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_ASS_CUUA_INVA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_ASS_CUUA_INVA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_ASS_CUUA_OTOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_ASS_CUUA_OTOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_ASS_CUUA_OTOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_ASS_CUUA_OTOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_ASS_CUUA_OTOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_ASS_CUUA_OTOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_ASS_CUUA_OTOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_ASS_CUUA_OTOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_ASS_CUUA_OTOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_ASS_CUUA_OTOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_ASS_CUUA_OTUA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_ASS_CUUA_OTUA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_ASS_CUUA_OTUA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_ASS_CUUA_OTUA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_ASS_CUUA_OTUA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_ASS_CUUA_OTUA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_ASS_CUUA_OTUA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_ASS_CUUA_OTUA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_ASS_CUUA_OTUA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_ASS_CUUA_OTUA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CGMB_ENABLED

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

#define CGMB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLDA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_EXP_CLDA_INVA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_EXP_CLDA_INVA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_EXP_CLDA_INVA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_EXP_CLDA_INVA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_EXP_CLDA_INVA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_EXP_CLDA_INVA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_EXP_CLDA_INVA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_EXP_CLDA_INVA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_EXP_CLDA_INVA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_EXP_CLDA_INVA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_EXP_CLDA_OTOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_EXP_CLDA_OTOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_EXP_CLDA_OTOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_EXP_CLDA_OTOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_EXP_CLDA_OTOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_EXP_CLDA_OTOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_EXP_CLDA_OTOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_EXP_CLDA_OTOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_EXP_CLDA_OTOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_EXP_CLDA_OTOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_EXP_CLDA_OTUA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_EXP_CLDA_OTUA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_EXP_CLDA_OTUA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_EXP_CLDA_OTUA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_EXP_CLDA_OTUA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_EXP_CLDA_OTUA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_EXP_CLDA_OTUA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_EXP_CLDA_OTUA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_EXP_CLDA_OTUA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_EXP_CLDA_OTUA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLDA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUDA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_EXP_CUDA_INVA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_EXP_CUDA_INVA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_EXP_CUDA_INVA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_EXP_CUDA_INVA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_EXP_CUDA_INVA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_EXP_CUDA_INVA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_EXP_CUDA_INVA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_EXP_CUDA_INVA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_EXP_CUDA_INVA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_EXP_CUDA_INVA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_EXP_CUDA_OTOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_EXP_CUDA_OTOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_EXP_CUDA_OTOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_EXP_CUDA_OTOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_EXP_CUDA_OTOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_EXP_CUDA_OTOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_EXP_CUDA_OTOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_EXP_CUDA_OTOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_EXP_CUDA_OTOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_EXP_CUDA_OTOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_EXP_CUDA_OTUA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_EXP_CUDA_OTUA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_EXP_CUDA_OTUA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_EXP_CUDA_OTUA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_EXP_CUDA_OTUA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_EXP_CUDA_OTUA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_EXP_CUDA_OTUA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_EXP_CUDA_OTUA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_EXP_CUDA_OTUA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_EXP_CUDA_OTUA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUDA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_EXP_CLUA_INVA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_EXP_CLUA_INVA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_EXP_CLUA_INVA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_EXP_CLUA_INVA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_EXP_CLUA_INVA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_EXP_CLUA_INVA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_EXP_CLUA_INVA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_EXP_CLUA_INVA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_EXP_CLUA_INVA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_EXP_CLUA_INVA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_EXP_CLUA_OTOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_EXP_CLUA_OTOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_EXP_CLUA_OTOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_EXP_CLUA_OTOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_EXP_CLUA_OTOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_EXP_CLUA_OTOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_EXP_CLUA_OTOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_EXP_CLUA_OTOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_EXP_CLUA_OTOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_EXP_CLUA_OTOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_EXP_CLUA_OTUA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_EXP_CLUA_OTUA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_EXP_CLUA_OTUA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_EXP_CLUA_OTUA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_EXP_CLUA_OTUA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_EXP_CLUA_OTUA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_EXP_CLUA_OTUA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_EXP_CLUA_OTUA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_EXP_CLUA_OTUA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_EXP_CLUA_OTUA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_EXP_CUUA_INVA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_EXP_CUUA_INVA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_EXP_CUUA_INVA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_EXP_CUUA_INVA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_EXP_CUUA_INVA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_EXP_CUUA_INVA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_EXP_CUUA_INVA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_EXP_CUUA_INVA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_EXP_CUUA_INVA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_EXP_CUUA_INVA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_EXP_CUUA_OTOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_EXP_CUUA_OTOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_EXP_CUUA_OTOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_EXP_CUUA_OTOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_EXP_CUUA_OTOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_EXP_CUUA_OTOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_EXP_CUUA_OTOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_EXP_CUUA_OTOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_EXP_CUUA_OTOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_EXP_CUUA_OTOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsv_EXP_CUUA_OTUA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsv_EXP_CUUA_OTUA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsv_EXP_CUUA_OTUA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsv_EXP_CUUA_OTUA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsv_EXP_CUUA_OTUA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsv_EXP_CUUA_OTUA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsv_EXP_CUUA_OTUA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsv_EXP_CUUA_OTUA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsv_EXP_CUUA_OTUA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsv_EXP_CUUA_OTUA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CGMB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef EXP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef trsv_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatMulTri_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatMulTri_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED
#define trmm_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ASS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CGMB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLDA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CLDA_ONOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CLDA_ONOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CLDA_ONOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CLDA_ONOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CLDA_ONOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CLDA_ONOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CLDA_ONOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CLDA_ONOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CLDA_ONOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CLDA_ONOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CLDA_OTSA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CLDA_OTSA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CLDA_OTSA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CLDA_OTSA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CLDA_OTSA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CLDA_OTSA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CLDA_OTSA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CLDA_OTSA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CLDA_OTSA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CLDA_OTSA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CLDA_OTHA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CLDA_OTHA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CLDA_OTHA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CLDA_OTHA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CLDA_OTHA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CLDA_OTHA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CLDA_OTHA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CLDA_OTHA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CLDA_OTHA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CLDA_OTHA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLDA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUDA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CUDA_ONOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CUDA_ONOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CUDA_ONOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CUDA_ONOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CUDA_ONOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CUDA_ONOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CUDA_ONOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CUDA_ONOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CUDA_ONOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CUDA_ONOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CUDA_OTSA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CUDA_OTSA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CUDA_OTSA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CUDA_OTSA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CUDA_OTSA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CUDA_OTSA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CUDA_OTSA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CUDA_OTSA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CUDA_OTSA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CUDA_OTSA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CUDA_OTHA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CUDA_OTHA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CUDA_OTHA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CUDA_OTHA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CUDA_OTHA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CUDA_OTHA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CUDA_OTHA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CUDA_OTHA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CUDA_OTHA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CUDA_OTHA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUDA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CLUA_ONOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CLUA_ONOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CLUA_ONOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CLUA_ONOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CLUA_ONOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CLUA_ONOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CLUA_ONOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CLUA_ONOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CLUA_ONOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CLUA_ONOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CLUA_OTSA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CLUA_OTSA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CLUA_OTSA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CLUA_OTSA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CLUA_OTSA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CLUA_OTSA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CLUA_OTSA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CLUA_OTSA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CLUA_OTSA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CLUA_OTSA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CLUA_OTHA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CLUA_OTHA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CLUA_OTHA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CLUA_OTHA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CLUA_OTHA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CLUA_OTHA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CLUA_OTHA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CLUA_OTHA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CLUA_OTHA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CLUA_OTHA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CUUA_ONOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CUUA_ONOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CUUA_ONOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CUUA_ONOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CUUA_ONOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CUUA_ONOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CUUA_ONOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CUUA_ONOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CUUA_ONOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CUUA_ONOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CUUA_OTSA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CUUA_OTSA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CUUA_OTSA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CUUA_OTSA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CUUA_OTSA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CUUA_OTSA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CUUA_OTSA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CUUA_OTSA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CUUA_OTSA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CUUA_OTSA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CUUA_OTHA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CUUA_OTHA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CUUA_OTHA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CUUA_OTHA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CUUA_OTHA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CUUA_OTHA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CUUA_OTHA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CUUA_OTHA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CUUA_OTHA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CUUA_OTHA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CGMB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CGMA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLDB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTSB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTSB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTSB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTSB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTSB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTSB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTSB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTSB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTSB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTSB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTHB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTHB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTHB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTHB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTHB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTHB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTHB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTHB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTHB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLDB_OTHB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLDB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUDB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTSB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTSB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTSB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTSB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTSB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTSB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTSB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTSB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTSB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTSB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTHB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTHB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTHB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTHB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTHB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTHB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTHB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTHB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTHB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUDB_OTHB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUDB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTSB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTSB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTSB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTSB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTSB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTSB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTSB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTSB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTSB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTSB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTHB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTHB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTHB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTHB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTHB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTHB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTHB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTHB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTHB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CLUB_OTHB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTSB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTSB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTSB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTSB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTSB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTSB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTSB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTSB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTSB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTSB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTHB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTHB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTHB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTHB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTHB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTHB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTHB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTHB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTHB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_ASS_CGMA_ONOA_CUUB_OTHB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CGMA_ENABLED

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

#define CGMB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLDA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CLDA_ONOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CLDA_ONOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CLDA_ONOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CLDA_ONOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CLDA_ONOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CLDA_ONOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CLDA_ONOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CLDA_ONOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CLDA_ONOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CLDA_ONOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CLDA_OTSA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CLDA_OTSA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CLDA_OTSA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CLDA_OTSA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CLDA_OTSA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CLDA_OTSA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CLDA_OTSA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CLDA_OTSA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CLDA_OTSA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CLDA_OTSA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CLDA_OTHA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CLDA_OTHA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CLDA_OTHA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CLDA_OTHA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CLDA_OTHA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CLDA_OTHA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CLDA_OTHA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CLDA_OTHA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CLDA_OTHA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CLDA_OTHA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLDA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUDA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CUDA_ONOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CUDA_ONOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CUDA_ONOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CUDA_ONOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CUDA_ONOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CUDA_ONOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CUDA_ONOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CUDA_ONOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CUDA_ONOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CUDA_ONOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CUDA_OTSA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CUDA_OTSA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CUDA_OTSA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CUDA_OTSA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CUDA_OTSA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CUDA_OTSA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CUDA_OTSA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CUDA_OTSA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CUDA_OTSA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CUDA_OTSA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CUDA_OTHA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CUDA_OTHA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CUDA_OTHA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CUDA_OTHA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CUDA_OTHA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CUDA_OTHA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CUDA_OTHA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CUDA_OTHA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CUDA_OTHA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CUDA_OTHA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUDA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CLUA_ONOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CLUA_ONOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CLUA_ONOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CLUA_ONOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CLUA_ONOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CLUA_ONOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CLUA_ONOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CLUA_ONOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CLUA_ONOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CLUA_ONOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CLUA_OTSA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CLUA_OTSA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CLUA_OTSA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CLUA_OTSA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CLUA_OTSA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CLUA_OTSA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CLUA_OTSA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CLUA_OTSA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CLUA_OTSA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CLUA_OTSA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CLUA_OTHA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CLUA_OTHA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CLUA_OTHA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CLUA_OTHA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CLUA_OTHA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CLUA_OTHA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CLUA_OTHA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CLUA_OTHA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CLUA_OTHA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CLUA_OTHA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CUUA_ONOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CUUA_ONOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CUUA_ONOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CUUA_ONOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CUUA_ONOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CUUA_ONOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CUUA_ONOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CUUA_ONOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CUUA_ONOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CUUA_ONOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CUUA_OTSA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CUUA_OTSA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CUUA_OTSA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CUUA_OTSA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CUUA_OTSA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CUUA_OTSA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CUUA_OTSA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CUUA_OTSA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CUUA_OTSA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CUUA_OTSA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CUUA_OTHA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CUUA_OTHA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CUUA_OTHA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CUUA_OTHA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CUUA_OTHA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CUUA_OTHA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CUUA_OTHA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CUUA_OTHA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CUUA_OTHA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CUUA_OTHA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CGMB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CGMA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLDB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTSB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTSB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTSB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTSB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTSB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTSB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTSB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTSB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTSB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTSB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTHB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTHB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTHB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTHB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTHB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTHB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTHB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTHB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTHB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLDB_OTHB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLDB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUDB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTSB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTSB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTSB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTSB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTSB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTSB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTSB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTSB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTSB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTSB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTHB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTHB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTHB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTHB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTHB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTHB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTHB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTHB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTHB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUDB_OTHB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUDB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTSB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTSB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTSB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTSB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTSB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTSB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTSB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTSB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTSB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTSB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTHB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTHB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTHB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTHB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTHB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTHB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTHB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTHB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTHB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CLUB_OTHB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTSB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTSB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTSB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTSB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTSB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTSB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTSB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTSB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTSB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTSB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTSB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTSB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTHB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTHB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTHB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTHB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTHB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTHB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTHB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTHB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTHB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTHB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trmm_EXP_CGMA_ONOA_CUUB_OTHB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTHB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CGMA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef EXP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef trmm_ENABLED
#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatMulTri_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatMulTri_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define trsm_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ASS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CGMB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLDA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CLDA_INVA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CLDA_INVA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CLDA_INVA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CLDA_INVA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CLDA_INVA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CLDA_INVA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CLDA_INVA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CLDA_INVA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CLDA_INVA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CLDA_INVA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CLDA_OTOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CLDA_OTOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CLDA_OTOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CLDA_OTOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CLDA_OTOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CLDA_OTOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CLDA_OTOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CLDA_OTOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CLDA_OTOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CLDA_OTOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CLDA_OTUA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CLDA_OTUA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CLDA_OTUA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CLDA_OTUA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CLDA_OTUA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CLDA_OTUA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CLDA_OTUA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CLDA_OTUA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CLDA_OTUA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CLDA_OTUA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLDA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUDA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CUDA_INVA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CUDA_INVA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CUDA_INVA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CUDA_INVA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CUDA_INVA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CUDA_INVA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CUDA_INVA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CUDA_INVA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CUDA_INVA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CUDA_INVA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CUDA_OTOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CUDA_OTOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CUDA_OTOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CUDA_OTOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CUDA_OTOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CUDA_OTOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CUDA_OTOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CUDA_OTOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CUDA_OTOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CUDA_OTOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CUDA_OTUA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CUDA_OTUA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CUDA_OTUA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CUDA_OTUA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CUDA_OTUA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CUDA_OTUA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CUDA_OTUA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CUDA_OTUA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CUDA_OTUA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CUDA_OTUA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUDA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CLUA_INVA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CLUA_INVA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CLUA_INVA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CLUA_INVA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CLUA_INVA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CLUA_INVA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CLUA_INVA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CLUA_INVA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CLUA_INVA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CLUA_INVA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CLUA_OTOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CLUA_OTOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CLUA_OTOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CLUA_OTOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CLUA_OTOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CLUA_OTOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CLUA_OTOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CLUA_OTOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CLUA_OTOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CLUA_OTOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CLUA_OTUA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CLUA_OTUA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CLUA_OTUA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CLUA_OTUA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CLUA_OTUA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CLUA_OTUA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CLUA_OTUA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CLUA_OTUA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CLUA_OTUA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CLUA_OTUA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CUUA_INVA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CUUA_INVA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CUUA_INVA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CUUA_INVA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CUUA_INVA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CUUA_INVA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CUUA_INVA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CUUA_INVA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CUUA_INVA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CUUA_INVA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CUUA_OTOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CUUA_OTOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CUUA_OTOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CUUA_OTOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CUUA_OTOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CUUA_OTOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CUUA_OTOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CUUA_OTOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CUUA_OTOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CUUA_OTOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CUUA_OTUA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CUUA_OTUA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CUUA_OTUA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CUUA_OTUA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CUUA_OTUA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CUUA_OTUA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CUUA_OTUA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CUUA_OTUA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CUUA_OTUA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CUUA_OTUA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CGMB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CGMA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLDB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_INVB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_INVB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_INVB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_INVB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_INVB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_INVB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_INVB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_INVB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_INVB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_INVB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTUB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTUB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTUB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTUB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTUB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTUB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTUB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTUB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTUB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLDB_OTUB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLDB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUDB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_INVB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_INVB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_INVB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_INVB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_INVB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_INVB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_INVB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_INVB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_INVB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_INVB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTUB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTUB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTUB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTUB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTUB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTUB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTUB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTUB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTUB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUDB_OTUB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUDB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_INVB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_INVB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_INVB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_INVB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_INVB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_INVB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_INVB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_INVB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_INVB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_INVB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTUB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTUB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTUB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTUB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTUB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTUB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTUB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTUB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTUB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CLUB_OTUB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_INVB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_INVB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_INVB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_INVB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_INVB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_INVB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_INVB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_INVB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_INVB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_INVB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTUB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTUB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTUB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTUB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTUB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTUB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTUB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTUB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTUB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_ASS_CGMA_ONOA_CUUB_OTUB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CGMA_ENABLED

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

#define CGMB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLDA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CLDA_INVA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CLDA_INVA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CLDA_INVA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CLDA_INVA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CLDA_INVA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CLDA_INVA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CLDA_INVA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CLDA_INVA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CLDA_INVA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CLDA_INVA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CLDA_OTOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CLDA_OTOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CLDA_OTOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CLDA_OTOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CLDA_OTOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CLDA_OTOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CLDA_OTOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CLDA_OTOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CLDA_OTOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CLDA_OTOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CLDA_OTUA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CLDA_OTUA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CLDA_OTUA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CLDA_OTUA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CLDA_OTUA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CLDA_OTUA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CLDA_OTUA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CLDA_OTUA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CLDA_OTUA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CLDA_OTUA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLDA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUDA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CUDA_INVA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CUDA_INVA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CUDA_INVA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CUDA_INVA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CUDA_INVA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CUDA_INVA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CUDA_INVA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CUDA_INVA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CUDA_INVA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CUDA_INVA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CUDA_OTOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CUDA_OTOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CUDA_OTOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CUDA_OTOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CUDA_OTOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CUDA_OTOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CUDA_OTOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CUDA_OTOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CUDA_OTOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CUDA_OTOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CUDA_OTUA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CUDA_OTUA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CUDA_OTUA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CUDA_OTUA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CUDA_OTUA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CUDA_OTUA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CUDA_OTUA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CUDA_OTUA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CUDA_OTUA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CUDA_OTUA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUDA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CLUA_INVA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CLUA_INVA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CLUA_INVA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CLUA_INVA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CLUA_INVA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CLUA_INVA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CLUA_INVA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CLUA_INVA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CLUA_INVA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CLUA_INVA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CLUA_OTOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CLUA_OTOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CLUA_OTOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CLUA_OTOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CLUA_OTOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CLUA_OTOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CLUA_OTOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CLUA_OTOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CLUA_OTOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CLUA_OTOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CLUA_OTUA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CLUA_OTUA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CLUA_OTUA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CLUA_OTUA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CLUA_OTUA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CLUA_OTUA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CLUA_OTUA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CLUA_OTUA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CLUA_OTUA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CLUA_OTUA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CUUA_INVA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CUUA_INVA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CUUA_INVA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CUUA_INVA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CUUA_INVA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CUUA_INVA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CUUA_INVA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CUUA_INVA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CUUA_INVA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CUUA_INVA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CUUA_OTOA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CUUA_OTOA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CUUA_OTOA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CUUA_OTOA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CUUA_OTOA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CUUA_OTOA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CUUA_OTOA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CUUA_OTOA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CUUA_OTOA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CUUA_OTOA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CUUA_OTUA_CGMB_ONOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CUUA_OTUA_CGMB_ONOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CUUA_OTUA_CGMB_ONOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CUUA_OTUA_CGMB_ONOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CUUA_OTUA_CGMB_ONOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CUUA_OTUA_CGMB_ONOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CUUA_OTUA_CGMB_ONOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CUUA_OTUA_CGMB_ONOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CUUA_OTUA_CGMB_ONOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CUUA_OTUA_CGMB_ONOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUUA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CGMB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CGMA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONOA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLDB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_INVB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_INVB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_INVB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_INVB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_INVB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_INVB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_INVB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_INVB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_INVB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_INVB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTUB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTUB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTUB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTUB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTUB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTUB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTUB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTUB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTUB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLDB_OTUB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLDB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUDB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_INVB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_INVB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_INVB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_INVB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_INVB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_INVB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_INVB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_INVB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_INVB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_INVB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTUB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTUB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTUB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTUB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTUB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTUB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTUB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTUB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTUB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUDB_OTUB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUDB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_INVB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_INVB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_INVB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_INVB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_INVB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_INVB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_INVB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_INVB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_INVB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_INVB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTUB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTUB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTUB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTUB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTUB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTUB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTUB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTUB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTUB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CLUB_OTUB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define INVB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_INVB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_INVB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_INVB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_INVB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_INVB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_INVB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_INVB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_INVB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_INVB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_INVB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef INVB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTOB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTOB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTOB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTOB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTOB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTOB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTOB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTOB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTOB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTOB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTOB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTOB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTUB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTUB_CK5
        use pm_kind, only: CKG => CK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTUB_CK4
        use pm_kind, only: CKG => CK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTUB_CK3
        use pm_kind, only: CKG => CK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTUB_CK2
        use pm_kind, only: CKG => CK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTUB_CK1
        use pm_kind, only: CKG => CK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTUB_RK5
        use pm_kind, only: RKG => RK5
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTUB_RK4
        use pm_kind, only: RKG => RK4
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTUB_RK3
        use pm_kind, only: RKG => RK3
#include "pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTUB_RK2
        use pm_kind, only: RKG => RK2
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure trsm_EXP_CGMA_ONOA_CUUB_OTUB_RK1
        use pm_kind, only: RKG => RK1
#define DISPATCH_ENABLED 1
#include "pm_matrixMulTri@routines.inc.F90"
#undef  DISPATCH_ENABLED
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUUB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONOA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CGMA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef EXP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef trsm_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatMulTri_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines