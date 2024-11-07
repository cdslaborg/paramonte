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

!>  \brief This file contains the implementations of the tests of module [pm_complexCompareAll](@ref pm_complexCompareAll).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_complexCompareAll) routines

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isallless_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isallless_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_isallless_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_isallless_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_isallless_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_isallless_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_isallless_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#undef isallless_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isallless_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isallleq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isallleq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_isallleq_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_isallleq_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_isallleq_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_isallleq_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_isallleq_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#undef isallleq_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isallleq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isallneq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isallneq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_isallneq_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_isallneq_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_isallneq_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_isallneq_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_isallneq_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#undef isallneq_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isallneq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isallmeq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isallmeq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_isallmeq_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_isallmeq_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_isallmeq_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_isallmeq_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_isallmeq_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#undef isallmeq_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isallmeq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isallmore_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isallmore_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_isallmore_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_isallmore_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_isallmore_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_isallmore_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_isallmore_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#undef isallmore_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isallmore_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE
