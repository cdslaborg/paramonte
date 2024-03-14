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

!>  \brief This file contains the implementations of the tests of module [val2complex_pmod](@ref pm_val2complex).
!>
!>  \fintest
!>
!>  \author
!>  \AmirShahmoradi

submodule (test_pm_val2complex) routines

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED
    module procedure test_getComplex128_LK_1
        use pm_val2complex, only: getComplex => getComplex128
        use pm_kind, only: IK, CK => CK3
#define test_getComplex128_LK_ENABLED 1
#include "test_pm_val2complex@routines.inc.F90"
#undef test_getComplex128_LK_ENABLED
    end procedure

    module procedure test_getComplex128_SK_1
        use pm_val2complex, only: getComplex => getComplex128
        use pm_kind, only: IK, CK => CK3
#define test_getComplex128_SK_ENABLED 1
#include "test_pm_val2complex@routines.inc.F90"
#undef test_getComplex128_SK_ENABLED
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getComplex64_LK_1
        use pm_val2complex, only: getComplex => getComplex64
        use pm_kind, only: IK, CK => CK2
#define test_getComplex64_LK_ENABLED 1
#include "test_pm_val2complex@routines.inc.F90"
#undef test_getComplex64_LK_ENABLED
    end procedure

    module procedure test_getComplex64_SK_1
        use pm_val2complex, only: getComplex => getComplex64
        use pm_kind, only: IK, CK => CK2
#define test_getComplex64_SK_ENABLED 1
#include "test_pm_val2complex@routines.inc.F90"
#undef test_getComplex64_SK_ENABLED
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getComplex32_LK_1
        use pm_val2complex, only: getComplex => getComplex32
        use pm_kind, only: IK, CK => CK1
#define test_getComplex32_LK_ENABLED 1
#include "test_pm_val2complex@routines.inc.F90"
#undef test_getComplex32_LK_ENABLED
    end procedure

    module procedure test_getComplex32_SK_1
        use pm_val2complex, only: getComplex => getComplex32
        use pm_kind, only: IK, CK => CK1
#define test_getComplex32_SK_ENABLED 1
#include "test_pm_val2complex@routines.inc.F90"
#undef test_getComplex32_SK_ENABLED
    end procedure
#endif

    module procedure test_getComplex_LK_1
        use pm_kind, only: IK, CK
        use pm_val2complex, only: getComplex
#define test_getComplex_LK_ENABLED 1
#include "test_pm_val2complex@routines.inc.F90"
#undef test_getComplex_LK_ENABLED
    end procedure

    module procedure test_getComplex_SK_1
        use pm_kind, only: IK, CK
        use pm_val2complex, only: getComplex
#define test_getComplex_SK_ENABLED 1
#include "test_pm_val2complex@routines.inc.F90"
#undef test_getComplex_SK_ENABLED
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE
