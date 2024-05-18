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
!>  This file contains procedure implementations of [val2complex_pmod](@ref pm_val2complex).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_val2complex) routines ! LCOV_EXCL_LINE

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getComplex_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getComplexDef_LK5
        use pm_kind, only: CKC => CK, LKC => LK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getComplexDef_LK4
        use pm_kind, only: CKC => CK, LKC => LK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getComplexDef_LK3
        use pm_kind, only: CKC => CK, LKC => LK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getComplexDef_LK2
        use pm_kind, only: CKC => CK, LKC => LK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getComplexDef_LK1
        use pm_kind, only: CKC => CK, LKC => LK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getComplexDef_SK5
        use pm_kind, only: CKC => CK, SKC => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getComplexDef_SK4
        use pm_kind, only: CKC => CK, SKC => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getComplexDef_SK3
        use pm_kind, only: CKC => CK, SKC => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getComplexDef_SK2
        use pm_kind, only: CKC => CK, SKC => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getComplexDef_SK1
        use pm_kind, only: CKC => CK, SKC => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getComplex_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setComplex_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED && LK5_ENABLED
    module procedure setComplexDef_CK5_LK5
        use pm_kind, only: CKC => CK5, LKC => LK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && LK4_ENABLED
    module procedure setComplexDef_CK5_LK4
        use pm_kind, only: CKC => CK5, LKC => LK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && LK3_ENABLED
    module procedure setComplexDef_CK5_LK3
        use pm_kind, only: CKC => CK5, LKC => LK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && LK2_ENABLED
    module procedure setComplexDef_CK5_LK2
        use pm_kind, only: CKC => CK5, LKC => LK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && LK1_ENABLED
    module procedure setComplexDef_CK5_LK1
        use pm_kind, only: CKC => CK5, LKC => LK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK4_ENABLED && LK5_ENABLED
    module procedure setComplexDef_CK4_LK5
        use pm_kind, only: CKC => CK4, LKC => LK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && LK4_ENABLED
    module procedure setComplexDef_CK4_LK4
        use pm_kind, only: CKC => CK4, LKC => LK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && LK3_ENABLED
    module procedure setComplexDef_CK4_LK3
        use pm_kind, only: CKC => CK4, LKC => LK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && LK2_ENABLED
    module procedure setComplexDef_CK4_LK2
        use pm_kind, only: CKC => CK4, LKC => LK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && LK1_ENABLED
    module procedure setComplexDef_CK4_LK1
        use pm_kind, only: CKC => CK4, LKC => LK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED && LK5_ENABLED
    module procedure setComplexDef_CK3_LK5
        use pm_kind, only: CKC => CK3, LKC => LK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && LK4_ENABLED
    module procedure setComplexDef_CK3_LK4
        use pm_kind, only: CKC => CK3, LKC => LK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && LK3_ENABLED
    module procedure setComplexDef_CK3_LK3
        use pm_kind, only: CKC => CK3, LKC => LK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && LK2_ENABLED
    module procedure setComplexDef_CK3_LK2
        use pm_kind, only: CKC => CK3, LKC => LK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && LK1_ENABLED
    module procedure setComplexDef_CK3_LK1
        use pm_kind, only: CKC => CK3, LKC => LK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK2_ENABLED && LK5_ENABLED
    module procedure setComplexDef_CK2_LK5
        use pm_kind, only: CKC => CK2, LKC => LK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && LK4_ENABLED
    module procedure setComplexDef_CK2_LK4
        use pm_kind, only: CKC => CK2, LKC => LK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && LK3_ENABLED
    module procedure setComplexDef_CK2_LK3
        use pm_kind, only: CKC => CK2, LKC => LK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && LK2_ENABLED
    module procedure setComplexDef_CK2_LK2
        use pm_kind, only: CKC => CK2, LKC => LK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && LK1_ENABLED
    module procedure setComplexDef_CK2_LK1
        use pm_kind, only: CKC => CK2, LKC => LK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK1_ENABLED && LK5_ENABLED
    module procedure setComplexDef_CK1_LK5
        use pm_kind, only: CKC => CK1, LKC => LK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && LK4_ENABLED
    module procedure setComplexDef_CK1_LK4
        use pm_kind, only: CKC => CK1, LKC => LK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && LK3_ENABLED
    module procedure setComplexDef_CK1_LK3
        use pm_kind, only: CKC => CK1, LKC => LK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && LK2_ENABLED
    module procedure setComplexDef_CK1_LK2
        use pm_kind, only: CKC => CK1, LKC => LK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && LK1_ENABLED
    module procedure setComplexDef_CK1_LK1
        use pm_kind, only: CKC => CK1, LKC => LK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED && SK5_ENABLED
    module procedure setComplexDef_CK5_SK5
        use pm_kind, only: CKC => CK5, SKC => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && SK4_ENABLED
    module procedure setComplexDef_CK5_SK4
        use pm_kind, only: CKC => CK5, SKC => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && SK3_ENABLED
    module procedure setComplexDef_CK5_SK3
        use pm_kind, only: CKC => CK5, SKC => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && SK2_ENABLED
    module procedure setComplexDef_CK5_SK2
        use pm_kind, only: CKC => CK5, SKC => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && SK1_ENABLED
    module procedure setComplexDef_CK5_SK1
        use pm_kind, only: CKC => CK5, SKC => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK4_ENABLED && SK5_ENABLED
    module procedure setComplexDef_CK4_SK5
        use pm_kind, only: CKC => CK4, SKC => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && SK4_ENABLED
    module procedure setComplexDef_CK4_SK4
        use pm_kind, only: CKC => CK4, SKC => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && SK3_ENABLED
    module procedure setComplexDef_CK4_SK3
        use pm_kind, only: CKC => CK4, SKC => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && SK2_ENABLED
    module procedure setComplexDef_CK4_SK2
        use pm_kind, only: CKC => CK4, SKC => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && SK1_ENABLED
    module procedure setComplexDef_CK4_SK1
        use pm_kind, only: CKC => CK4, SKC => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED && SK5_ENABLED
    module procedure setComplexDef_CK3_SK5
        use pm_kind, only: CKC => CK3, SKC => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && SK4_ENABLED
    module procedure setComplexDef_CK3_SK4
        use pm_kind, only: CKC => CK3, SKC => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && SK3_ENABLED
    module procedure setComplexDef_CK3_SK3
        use pm_kind, only: CKC => CK3, SKC => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && SK2_ENABLED
    module procedure setComplexDef_CK3_SK2
        use pm_kind, only: CKC => CK3, SKC => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && SK1_ENABLED
    module procedure setComplexDef_CK3_SK1
        use pm_kind, only: CKC => CK3, SKC => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK2_ENABLED && SK5_ENABLED
    module procedure setComplexDef_CK2_SK5
        use pm_kind, only: CKC => CK2, SKC => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && SK4_ENABLED
    module procedure setComplexDef_CK2_SK4
        use pm_kind, only: CKC => CK2, SKC => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && SK3_ENABLED
    module procedure setComplexDef_CK2_SK3
        use pm_kind, only: CKC => CK2, SKC => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && SK2_ENABLED
    module procedure setComplexDef_CK2_SK2
        use pm_kind, only: CKC => CK2, SKC => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && SK1_ENABLED
    module procedure setComplexDef_CK2_SK1
        use pm_kind, only: CKC => CK2, SKC => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK1_ENABLED && SK5_ENABLED
    module procedure setComplexDef_CK1_SK5
        use pm_kind, only: CKC => CK1, SKC => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && SK4_ENABLED
    module procedure setComplexDef_CK1_SK4
        use pm_kind, only: CKC => CK1, SKC => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && SK3_ENABLED
    module procedure setComplexDef_CK1_SK3
        use pm_kind, only: CKC => CK1, SKC => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && SK2_ENABLED
    module procedure setComplexDef_CK1_SK2
        use pm_kind, only: CKC => CK1, SKC => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && SK1_ENABLED
    module procedure setComplexDef_CK1_SK1
        use pm_kind, only: CKC => CK1, SKC => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Err_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED && SK5_ENABLED
    module procedure setComplexIO_CK5_SK5
        use pm_kind, only: CKC => CK5, SKC => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && SK4_ENABLED
    module procedure setComplexIO_CK5_SK4
        use pm_kind, only: CKC => CK5, SKC => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && SK3_ENABLED
    module procedure setComplexIO_CK5_SK3
        use pm_kind, only: CKC => CK5, SKC => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && SK2_ENABLED
    module procedure setComplexIO_CK5_SK2
        use pm_kind, only: CKC => CK5, SKC => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && SK1_ENABLED
    module procedure setComplexIO_CK5_SK1
        use pm_kind, only: CKC => CK5, SKC => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK4_ENABLED && SK5_ENABLED
    module procedure setComplexIO_CK4_SK5
        use pm_kind, only: CKC => CK4, SKC => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && SK4_ENABLED
    module procedure setComplexIO_CK4_SK4
        use pm_kind, only: CKC => CK4, SKC => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && SK3_ENABLED
    module procedure setComplexIO_CK4_SK3
        use pm_kind, only: CKC => CK4, SKC => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && SK2_ENABLED
    module procedure setComplexIO_CK4_SK2
        use pm_kind, only: CKC => CK4, SKC => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && SK1_ENABLED
    module procedure setComplexIO_CK4_SK1
        use pm_kind, only: CKC => CK4, SKC => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED && SK5_ENABLED
    module procedure setComplexIO_CK3_SK5
        use pm_kind, only: CKC => CK3, SKC => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && SK4_ENABLED
    module procedure setComplexIO_CK3_SK4
        use pm_kind, only: CKC => CK3, SKC => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && SK3_ENABLED
    module procedure setComplexIO_CK3_SK3
        use pm_kind, only: CKC => CK3, SKC => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && SK2_ENABLED
    module procedure setComplexIO_CK3_SK2
        use pm_kind, only: CKC => CK3, SKC => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && SK1_ENABLED
    module procedure setComplexIO_CK3_SK1
        use pm_kind, only: CKC => CK3, SKC => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK2_ENABLED && SK5_ENABLED
    module procedure setComplexIO_CK2_SK5
        use pm_kind, only: CKC => CK2, SKC => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && SK4_ENABLED
    module procedure setComplexIO_CK2_SK4
        use pm_kind, only: CKC => CK2, SKC => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && SK3_ENABLED
    module procedure setComplexIO_CK2_SK3
        use pm_kind, only: CKC => CK2, SKC => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && SK2_ENABLED
    module procedure setComplexIO_CK2_SK2
        use pm_kind, only: CKC => CK2, SKC => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && SK1_ENABLED
    module procedure setComplexIO_CK2_SK1
        use pm_kind, only: CKC => CK2, SKC => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK1_ENABLED && SK5_ENABLED
    module procedure setComplexIO_CK1_SK5
        use pm_kind, only: CKC => CK1, SKC => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && SK4_ENABLED
    module procedure setComplexIO_CK1_SK4
        use pm_kind, only: CKC => CK1, SKC => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && SK3_ENABLED
    module procedure setComplexIO_CK1_SK3
        use pm_kind, only: CKC => CK1, SKC => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && SK2_ENABLED
    module procedure setComplexIO_CK1_SK2
        use pm_kind, only: CKC => CK1, SKC => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && SK1_ENABLED
    module procedure setComplexIO_CK1_SK1
        use pm_kind, only: CKC => CK1, SKC => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Err_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setComplex_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
