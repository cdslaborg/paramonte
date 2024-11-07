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
!>  This file contains procedure implementations of [pm_complexCompareAll](@ref pm_complexCompareAll).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_complexCompareAll) routines ! LCOV_EXCL_LINE

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isallless_CK_ENABLED 1

#if CK5_ENABLED
    module procedure isallless_CK5
        use pm_kind, only: CKG => CK5
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isallless_CK4
        use pm_kind, only: CKG => CK4
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isallless_CK3
        use pm_kind, only: CKG => CK3
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isallless_CK2
        use pm_kind, only: CKG => CK2
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isallless_CK1
        use pm_kind, only: CKG => CK1
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#undef isallless_CK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isallleq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure isallleq_CK5
        use pm_kind, only: CKG => CK5
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isallleq_CK4
        use pm_kind, only: CKG => CK4
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isallleq_CK3
        use pm_kind, only: CKG => CK3
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isallleq_CK2
        use pm_kind, only: CKG => CK2
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isallleq_CK1
        use pm_kind, only: CKG => CK1
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#undef isallleq_CK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isallneq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure isallneq_CK5
        use pm_kind, only: CKG => CK5
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isallneq_CK4
        use pm_kind, only: CKG => CK4
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isallneq_CK3
        use pm_kind, only: CKG => CK3
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isallneq_CK2
        use pm_kind, only: CKG => CK2
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isallneq_CK1
        use pm_kind, only: CKG => CK1
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#undef isallneq_CK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define isalleq_CK_ENABLED 1
!
!#if CK3_ENABLED
!    module procedure isalleq_CK3
!        use pm_kind, only: CKG => CK3
!#include "pm_complexCompareAll@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK2_ENABLED
!    module procedure isalleq_CK2
!        use pm_kind, only: CKG => CK2
!#include "pm_complexCompareAll@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK1_ENABLED
!    module procedure isalleq_CK1
!        use pm_kind, only: CKG => CK1
!#include "pm_complexCompareAll@routines.inc.F90"
!    end procedure
!#endif
!
!#undef isalleq_CK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isallmeq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure isallmeq_CK5
        use pm_kind, only: CKG => CK5
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isallmeq_CK4
        use pm_kind, only: CKG => CK4
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isallmeq_CK3
        use pm_kind, only: CKG => CK3
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isallmeq_CK2
        use pm_kind, only: CKG => CK2
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isallmeq_CK1
        use pm_kind, only: CKG => CK1
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#undef isallmeq_CK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isallmore_CK_ENABLED 1

#if CK5_ENABLED
    module procedure isallmore_CK5
        use pm_kind, only: CKG => CK5
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isallmore_CK4
        use pm_kind, only: CKG => CK4
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isallmore_CK3
        use pm_kind, only: CKG => CK3
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isallmore_CK2
        use pm_kind, only: CKG => CK2
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isallmore_CK1
        use pm_kind, only: CKG => CK1
#include "pm_complexCompareAll@routines.inc.F90"
    end procedure
#endif

#undef isallmore_CK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
