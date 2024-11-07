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
!>  This file contains procedure implementations of [pm_complexCompareAny](@ref pm_complexCompareAny).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_complexCompareAny) routines ! LCOV_EXCL_LINE

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanyless_CK_ENABLED 1

#if CK5_ENABLED
    module procedure isanyless_CK5
        use pm_kind, only: CKG => CK5
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isanyless_CK4
        use pm_kind, only: CKG => CK4
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isanyless_CK3
        use pm_kind, only: CKG => CK3
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isanyless_CK2
        use pm_kind, only: CKG => CK2
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isanyless_CK1
        use pm_kind, only: CKG => CK1
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#undef isanyless_CK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanyleq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure isanyleq_CK5
        use pm_kind, only: CKG => CK5
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isanyleq_CK4
        use pm_kind, only: CKG => CK4
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isanyleq_CK3
        use pm_kind, only: CKG => CK3
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isanyleq_CK2
        use pm_kind, only: CKG => CK2
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isanyleq_CK1
        use pm_kind, only: CKG => CK1
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#undef isanyleq_CK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanyneq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure isanyneq_CK5
        use pm_kind, only: CKG => CK5
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isanyneq_CK4
        use pm_kind, only: CKG => CK4
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isanyneq_CK3
        use pm_kind, only: CKG => CK3
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isanyneq_CK2
        use pm_kind, only: CKG => CK2
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isanyneq_CK1
        use pm_kind, only: CKG => CK1
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#undef isanyneq_CK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanyeq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure isanyeq_CK5
        use pm_kind, only: CKG => CK5
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isanyeq_CK4
        use pm_kind, only: CKG => CK4
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isanyeq_CK3
        use pm_kind, only: CKG => CK3
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isanyeq_CK2
        use pm_kind, only: CKG => CK2
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isanyeq_CK1
        use pm_kind, only: CKG => CK1
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#undef isanyeq_CK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanymeq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure isanymeq_CK5
        use pm_kind, only: CKG => CK5
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isanymeq_CK4
        use pm_kind, only: CKG => CK4
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isanymeq_CK3
        use pm_kind, only: CKG => CK3
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isanymeq_CK2
        use pm_kind, only: CKG => CK2
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isanymeq_CK1
        use pm_kind, only: CKG => CK1
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#undef isanymeq_CK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isanymore_CK_ENABLED 1

#if CK5_ENABLED
    module procedure isanymore_CK5
        use pm_kind, only: CKG => CK5
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure isanymore_CK4
        use pm_kind, only: CKG => CK4
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure isanymore_CK3
        use pm_kind, only: CKG => CK3
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure isanymore_CK2
        use pm_kind, only: CKG => CK2
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure isanymore_CK1
        use pm_kind, only: CKG => CK1
#include "pm_complexCompareAny@routines.inc.F90"
    end procedure
#endif

#undef isanymore_CK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
