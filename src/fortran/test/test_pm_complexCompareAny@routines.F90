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

!>  \brief This file contains the implementations of the tests of module [pm_complexCompareAny](@ref pm_complexCompareAny).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_complexCompareAny) routines

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanyless_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanyless_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_isanyless_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_isanyless_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_isanyless_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_isanyless_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_isanyless_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#undef isanyless_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isanyless_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanyleq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanyleq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_isanyleq_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_isanyleq_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_isanyleq_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_isanyleq_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_isanyleq_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#undef isanyleq_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isanyleq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanyneq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanyneq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_isanyneq_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_isanyneq_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_isanyneq_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_isanyneq_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_isanyneq_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#undef isanyneq_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isanyneq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanyeq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanyeq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_isanyeq_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_isanyeq_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_isanyeq_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_isanyeq_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_isanyeq_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#undef isanyeq_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isanyeq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanymeq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanymeq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_isanymeq_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_isanymeq_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_isanymeq_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_isanymeq_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_isanymeq_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#undef isanymeq_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isanymeq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanymore_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanymore_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_isanymore_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_isanymore_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_isanymore_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_isanymore_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_isanymore_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#undef isanymore_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isanymore_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE
