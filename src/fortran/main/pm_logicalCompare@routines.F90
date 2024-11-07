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
!>  This file contains procedure implementations of [pm_logicalCompare](@ref pm_logicalCompare).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_logicalCompare) routines ! LCOV_EXCL_LINE

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isless_ENABLED 1

#if LK5_ENABLED
    module procedure isless_LK5
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure isless_LK4
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure isless_LK3
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure isless_LK2
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure isless_LK1
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#undef isless_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isleq_ENABLED 1

#if LK5_ENABLED
    module procedure isleq_LK5
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure isleq_LK4
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure isleq_LK3
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure isleq_LK2
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure isleq_LK1
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#undef isleq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define iseq_ENABLED 1

#if LK5_ENABLED
    module procedure iseq_LK5
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure iseq_LK4
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure iseq_LK3
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure iseq_LK2
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure iseq_LK1
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#undef iseq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isneq_ENABLED 1

#if LK5_ENABLED
    module procedure isneq_LK5
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure isneq_LK4
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure isneq_LK3
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure isneq_LK2
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure isneq_LK1
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#undef isneq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ismeq_ENABLED 1

#if LK5_ENABLED
    module procedure ismeq_LK5
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure ismeq_LK4
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure ismeq_LK3
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure ismeq_LK2
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure ismeq_LK1
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#undef ismeq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ismore_ENABLED 1

#if LK5_ENABLED
    module procedure ismore_LK5
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure ismore_LK4
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure ismore_LK3
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure ismore_LK2
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure ismore_LK1
#include "pm_logicalCompare@routines.inc.F90"
    end procedure
#endif

#undef ismore_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines