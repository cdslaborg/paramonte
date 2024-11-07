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

!>  \brief This file contains the implementations of the tests of module [pm_complexCompareLex](@ref pm_complexCompareLex).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

submodule (test_pm_complexCompareLex) routines

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define islexless_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define islexless_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_islexless_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_islexless_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_islexless_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_islexless_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_islexless_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#undef islexless_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef islexless_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define islexleq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define islexleq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_islexleq_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_islexleq_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_islexleq_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_islexleq_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_islexleq_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#undef islexleq_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef islexleq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define islexmeq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define islexmeq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_islexmeq_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_islexmeq_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_islexmeq_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_islexmeq_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_islexmeq_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#undef islexmeq_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef islexmeq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define islexmore_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define islexmore_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_islexmore_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_islexmore_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_islexmore_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_islexmore_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_islexmore_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#undef islexmore_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef islexmore_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE
