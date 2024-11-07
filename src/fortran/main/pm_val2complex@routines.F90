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
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

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
        use pm_kind, only: CKG => CK, LKG => LK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getComplexDef_LK4
        use pm_kind, only: CKG => CK, LKG => LK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getComplexDef_LK3
        use pm_kind, only: CKG => CK, LKG => LK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getComplexDef_LK2
        use pm_kind, only: CKG => CK, LKG => LK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getComplexDef_LK1
        use pm_kind, only: CKG => CK, LKG => LK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getComplexDef_SK5
        use pm_kind, only: CKG => CK, SKG => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getComplexDef_SK4
        use pm_kind, only: CKG => CK, SKG => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getComplexDef_SK3
        use pm_kind, only: CKG => CK, SKG => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getComplexDef_SK2
        use pm_kind, only: CKG => CK, SKG => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getComplexDef_SK1
        use pm_kind, only: CKG => CK, SKG => SK1
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
        use pm_kind, only: CKG => CK5, LKG => LK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && LK4_ENABLED
    module procedure setComplexDef_CK5_LK4
        use pm_kind, only: CKG => CK5, LKG => LK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && LK3_ENABLED
    module procedure setComplexDef_CK5_LK3
        use pm_kind, only: CKG => CK5, LKG => LK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && LK2_ENABLED
    module procedure setComplexDef_CK5_LK2
        use pm_kind, only: CKG => CK5, LKG => LK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && LK1_ENABLED
    module procedure setComplexDef_CK5_LK1
        use pm_kind, only: CKG => CK5, LKG => LK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK4_ENABLED && LK5_ENABLED
    module procedure setComplexDef_CK4_LK5
        use pm_kind, only: CKG => CK4, LKG => LK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && LK4_ENABLED
    module procedure setComplexDef_CK4_LK4
        use pm_kind, only: CKG => CK4, LKG => LK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && LK3_ENABLED
    module procedure setComplexDef_CK4_LK3
        use pm_kind, only: CKG => CK4, LKG => LK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && LK2_ENABLED
    module procedure setComplexDef_CK4_LK2
        use pm_kind, only: CKG => CK4, LKG => LK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && LK1_ENABLED
    module procedure setComplexDef_CK4_LK1
        use pm_kind, only: CKG => CK4, LKG => LK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED && LK5_ENABLED
    module procedure setComplexDef_CK3_LK5
        use pm_kind, only: CKG => CK3, LKG => LK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && LK4_ENABLED
    module procedure setComplexDef_CK3_LK4
        use pm_kind, only: CKG => CK3, LKG => LK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && LK3_ENABLED
    module procedure setComplexDef_CK3_LK3
        use pm_kind, only: CKG => CK3, LKG => LK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && LK2_ENABLED
    module procedure setComplexDef_CK3_LK2
        use pm_kind, only: CKG => CK3, LKG => LK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && LK1_ENABLED
    module procedure setComplexDef_CK3_LK1
        use pm_kind, only: CKG => CK3, LKG => LK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK2_ENABLED && LK5_ENABLED
    module procedure setComplexDef_CK2_LK5
        use pm_kind, only: CKG => CK2, LKG => LK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && LK4_ENABLED
    module procedure setComplexDef_CK2_LK4
        use pm_kind, only: CKG => CK2, LKG => LK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && LK3_ENABLED
    module procedure setComplexDef_CK2_LK3
        use pm_kind, only: CKG => CK2, LKG => LK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && LK2_ENABLED
    module procedure setComplexDef_CK2_LK2
        use pm_kind, only: CKG => CK2, LKG => LK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && LK1_ENABLED
    module procedure setComplexDef_CK2_LK1
        use pm_kind, only: CKG => CK2, LKG => LK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK1_ENABLED && LK5_ENABLED
    module procedure setComplexDef_CK1_LK5
        use pm_kind, only: CKG => CK1, LKG => LK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && LK4_ENABLED
    module procedure setComplexDef_CK1_LK4
        use pm_kind, only: CKG => CK1, LKG => LK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && LK3_ENABLED
    module procedure setComplexDef_CK1_LK3
        use pm_kind, only: CKG => CK1, LKG => LK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && LK2_ENABLED
    module procedure setComplexDef_CK1_LK2
        use pm_kind, only: CKG => CK1, LKG => LK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && LK1_ENABLED
    module procedure setComplexDef_CK1_LK1
        use pm_kind, only: CKG => CK1, LKG => LK1
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
        use pm_kind, only: CKG => CK5, SKG => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && SK4_ENABLED
    module procedure setComplexDef_CK5_SK4
        use pm_kind, only: CKG => CK5, SKG => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && SK3_ENABLED
    module procedure setComplexDef_CK5_SK3
        use pm_kind, only: CKG => CK5, SKG => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && SK2_ENABLED
    module procedure setComplexDef_CK5_SK2
        use pm_kind, only: CKG => CK5, SKG => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && SK1_ENABLED
    module procedure setComplexDef_CK5_SK1
        use pm_kind, only: CKG => CK5, SKG => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK4_ENABLED && SK5_ENABLED
    module procedure setComplexDef_CK4_SK5
        use pm_kind, only: CKG => CK4, SKG => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && SK4_ENABLED
    module procedure setComplexDef_CK4_SK4
        use pm_kind, only: CKG => CK4, SKG => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && SK3_ENABLED
    module procedure setComplexDef_CK4_SK3
        use pm_kind, only: CKG => CK4, SKG => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && SK2_ENABLED
    module procedure setComplexDef_CK4_SK2
        use pm_kind, only: CKG => CK4, SKG => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && SK1_ENABLED
    module procedure setComplexDef_CK4_SK1
        use pm_kind, only: CKG => CK4, SKG => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED && SK5_ENABLED
    module procedure setComplexDef_CK3_SK5
        use pm_kind, only: CKG => CK3, SKG => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && SK4_ENABLED
    module procedure setComplexDef_CK3_SK4
        use pm_kind, only: CKG => CK3, SKG => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && SK3_ENABLED
    module procedure setComplexDef_CK3_SK3
        use pm_kind, only: CKG => CK3, SKG => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && SK2_ENABLED
    module procedure setComplexDef_CK3_SK2
        use pm_kind, only: CKG => CK3, SKG => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && SK1_ENABLED
    module procedure setComplexDef_CK3_SK1
        use pm_kind, only: CKG => CK3, SKG => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK2_ENABLED && SK5_ENABLED
    module procedure setComplexDef_CK2_SK5
        use pm_kind, only: CKG => CK2, SKG => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && SK4_ENABLED
    module procedure setComplexDef_CK2_SK4
        use pm_kind, only: CKG => CK2, SKG => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && SK3_ENABLED
    module procedure setComplexDef_CK2_SK3
        use pm_kind, only: CKG => CK2, SKG => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && SK2_ENABLED
    module procedure setComplexDef_CK2_SK2
        use pm_kind, only: CKG => CK2, SKG => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && SK1_ENABLED
    module procedure setComplexDef_CK2_SK1
        use pm_kind, only: CKG => CK2, SKG => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK1_ENABLED && SK5_ENABLED
    module procedure setComplexDef_CK1_SK5
        use pm_kind, only: CKG => CK1, SKG => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && SK4_ENABLED
    module procedure setComplexDef_CK1_SK4
        use pm_kind, only: CKG => CK1, SKG => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && SK3_ENABLED
    module procedure setComplexDef_CK1_SK3
        use pm_kind, only: CKG => CK1, SKG => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && SK2_ENABLED
    module procedure setComplexDef_CK1_SK2
        use pm_kind, only: CKG => CK1, SKG => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && SK1_ENABLED
    module procedure setComplexDef_CK1_SK1
        use pm_kind, only: CKG => CK1, SKG => SK1
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
        use pm_kind, only: CKG => CK5, SKG => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && SK4_ENABLED
    module procedure setComplexIO_CK5_SK4
        use pm_kind, only: CKG => CK5, SKG => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && SK3_ENABLED
    module procedure setComplexIO_CK5_SK3
        use pm_kind, only: CKG => CK5, SKG => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && SK2_ENABLED
    module procedure setComplexIO_CK5_SK2
        use pm_kind, only: CKG => CK5, SKG => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED && SK1_ENABLED
    module procedure setComplexIO_CK5_SK1
        use pm_kind, only: CKG => CK5, SKG => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK4_ENABLED && SK5_ENABLED
    module procedure setComplexIO_CK4_SK5
        use pm_kind, only: CKG => CK4, SKG => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && SK4_ENABLED
    module procedure setComplexIO_CK4_SK4
        use pm_kind, only: CKG => CK4, SKG => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && SK3_ENABLED
    module procedure setComplexIO_CK4_SK3
        use pm_kind, only: CKG => CK4, SKG => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && SK2_ENABLED
    module procedure setComplexIO_CK4_SK2
        use pm_kind, only: CKG => CK4, SKG => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED && SK1_ENABLED
    module procedure setComplexIO_CK4_SK1
        use pm_kind, only: CKG => CK4, SKG => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED && SK5_ENABLED
    module procedure setComplexIO_CK3_SK5
        use pm_kind, only: CKG => CK3, SKG => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && SK4_ENABLED
    module procedure setComplexIO_CK3_SK4
        use pm_kind, only: CKG => CK3, SKG => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && SK3_ENABLED
    module procedure setComplexIO_CK3_SK3
        use pm_kind, only: CKG => CK3, SKG => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && SK2_ENABLED
    module procedure setComplexIO_CK3_SK2
        use pm_kind, only: CKG => CK3, SKG => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED && SK1_ENABLED
    module procedure setComplexIO_CK3_SK1
        use pm_kind, only: CKG => CK3, SKG => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK2_ENABLED && SK5_ENABLED
    module procedure setComplexIO_CK2_SK5
        use pm_kind, only: CKG => CK2, SKG => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && SK4_ENABLED
    module procedure setComplexIO_CK2_SK4
        use pm_kind, only: CKG => CK2, SKG => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && SK3_ENABLED
    module procedure setComplexIO_CK2_SK3
        use pm_kind, only: CKG => CK2, SKG => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && SK2_ENABLED
    module procedure setComplexIO_CK2_SK2
        use pm_kind, only: CKG => CK2, SKG => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED && SK1_ENABLED
    module procedure setComplexIO_CK2_SK1
        use pm_kind, only: CKG => CK2, SKG => SK1
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK1_ENABLED && SK5_ENABLED
    module procedure setComplexIO_CK1_SK5
        use pm_kind, only: CKG => CK1, SKG => SK5
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && SK4_ENABLED
    module procedure setComplexIO_CK1_SK4
        use pm_kind, only: CKG => CK1, SKG => SK4
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && SK3_ENABLED
    module procedure setComplexIO_CK1_SK3
        use pm_kind, only: CKG => CK1, SKG => SK3
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && SK2_ENABLED
    module procedure setComplexIO_CK1_SK2
        use pm_kind, only: CKG => CK1, SKG => SK2
#include "pm_val2complex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED && SK1_ENABLED
    module procedure setComplexIO_CK1_SK1
        use pm_kind, only: CKG => CK1, SKG => SK1
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
