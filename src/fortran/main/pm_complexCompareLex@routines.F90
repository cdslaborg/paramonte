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
!>  This file contains procedure implementations of [pm_complexCompareLex](@ref pm_complexCompareLex).
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

submodule (pm_complexCompareLex) routines ! LCOV_EXCL_LINE

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define islexless_CK_ENABLED 1

#if CK5_ENABLED
    module procedure islexless_CK5
        use pm_kind, only: CKG => CK5
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure islexless_CK4
        use pm_kind, only: CKG => CK4
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure islexless_CK3
        use pm_kind, only: CKG => CK3
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure islexless_CK2
        use pm_kind, only: CKG => CK2
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure islexless_CK1
        use pm_kind, only: CKG => CK1
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#undef islexless_CK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define islexleq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure islexleq_CK5
        use pm_kind, only: CKG => CK5
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure islexleq_CK4
        use pm_kind, only: CKG => CK4
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure islexleq_CK3
        use pm_kind, only: CKG => CK3
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure islexleq_CK2
        use pm_kind, only: CKG => CK2
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure islexleq_CK1
        use pm_kind, only: CKG => CK1
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#undef islexleq_CK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define islexneq_CK_ENABLED 1
!
!#if CK3_ENABLED
!    module procedure islexneq_CK3
!        use pm_kind, only: CKG => CK3
!#include "pm_complexCompareLex@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK2_ENABLED
!    module procedure islexneq_CK2
!        use pm_kind, only: CKG => CK2
!#include "pm_complexCompareLex@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK1_ENABLED
!    module procedure islexneq_CK1
!        use pm_kind, only: CKG => CK1
!#include "pm_complexCompareLex@routines.inc.F90"
!    end procedure
!#endif
!
!#undef islexneq_CK_ENABLED
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define isanyeq_CK_ENABLED 1
!
!#if CK3_ENABLED
!    module procedure isanyeq_CK3
!        use pm_kind, only: CKG => CK3
!#include "pm_complexCompareLex@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK2_ENABLED
!    module procedure isanyeq_CK2
!        use pm_kind, only: CKG => CK2
!#include "pm_complexCompareLex@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK1_ENABLED
!    module procedure isanyeq_CK1
!        use pm_kind, only: CKG => CK1
!#include "pm_complexCompareLex@routines.inc.F90"
!    end procedure
!#endif
!
!#undef isanyeq_CK_ENABLED
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define islexmeq_CK_ENABLED 1

#if CK5_ENABLED
    module procedure islexmeq_CK5
        use pm_kind, only: CKG => CK5
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure islexmeq_CK4
        use pm_kind, only: CKG => CK4
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure islexmeq_CK3
        use pm_kind, only: CKG => CK3
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure islexmeq_CK2
        use pm_kind, only: CKG => CK2
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure islexmeq_CK1
        use pm_kind, only: CKG => CK1
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#undef islexmeq_CK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define islexmore_CK_ENABLED 1

#if CK5_ENABLED
    module procedure islexmore_CK5
        use pm_kind, only: CKG => CK5
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure islexmore_CK4
        use pm_kind, only: CKG => CK4
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure islexmore_CK3
        use pm_kind, only: CKG => CK3
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure islexmore_CK2
        use pm_kind, only: CKG => CK2
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure islexmore_CK1
        use pm_kind, only: CKG => CK1
#include "pm_complexCompareLex@routines.inc.F90"
    end procedure
#endif

#undef islexmore_CK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
