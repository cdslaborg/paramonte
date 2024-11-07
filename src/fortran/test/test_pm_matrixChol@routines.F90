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
!>  This file contains the implementations of the tests of module [pm_matrixChol](@ref pm_matrixChol).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_matrixChol) routines

    use pm_err, only: getFine
    use pm_distCov, only: getCovRand
    use pm_arrayFill, only: getFilled
    use pm_distUnif, only: getUnifRand
    use pm_matrixCopy, only: getMatCopy, setMatCopy
    use pm_matrixInit, only: setMatInit, upp, dia, low
    use pm_matrixInit, only: getMatInit, uppLowDia
    use pm_complexCompareAll, only: operator(<)
    use pm_matrixSubset, only: subset_type
    use pm_arrayResize, only: setResized
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setChoLow_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setChoLow_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setChoLow_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setChoLow_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setChoLow_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setChoLow_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#undef  RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  setChoLow_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatChol_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getMatChol_CK5
        use pm_kind, only: TKG => CK5
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getMatChol_CK4
        use pm_kind, only: TKG => CK4
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getMatChol_CK3
        use pm_kind, only: TKG => CK3
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getMatChol_CK2
        use pm_kind, only: TKG => CK2
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getMatChol_CK1
        use pm_kind, only: TKG => CK1
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#undef  CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getMatChol_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getMatChol_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getMatChol_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getMatChol_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getMatChol_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#undef  RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  getMatChol_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatChol_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setMatChol_CK5
        use pm_kind, only: TKG => CK5
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setMatChol_CK4
        use pm_kind, only: TKG => CK4
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setMatChol_CK3
        use pm_kind, only: TKG => CK3
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setMatChol_CK2
        use pm_kind, only: TKG => CK2
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setMatChol_CK1
        use pm_kind, only: TKG => CK1
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#undef  CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setMatChol_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setMatChol_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setMatChol_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setMatChol_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setMatChol_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_matrixChol@routines.inc.F90"
    end procedure
#endif

#undef  RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  setMatChol_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE