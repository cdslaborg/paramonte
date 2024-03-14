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

!>  \brief This file contains the implementations of the tests of module [pm_matrixMulTri](@ref pm_matrixMulTri).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_matrixMulTri) routines

    use pm_kind, only: RKB
    use pm_option, only: getOption
    use pm_io, only: display_type
    use pm_mathCompare, only: isClose
    use pm_arrayChoice, only: getChoice
    use pm_distUnif, only: getUnifRand, setUnifRand

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatMulTri_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setMatMulTri_CK5
        use pm_kind, only: TKC => CK5
#include "test_pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setMatMulTri_CK4
        use pm_kind, only: TKC => CK4
#include "test_pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setMatMulTri_CK3
        use pm_kind, only: TKC => CK3
#include "test_pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setMatMulTri_CK2
        use pm_kind, only: TKC => CK2
#include "test_pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setMatMulTri_CK1
        use pm_kind, only: TKC => CK1
#include "test_pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setMatMulTri_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setMatMulTri_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setMatMulTri_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setMatMulTri_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setMatMulTri_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_matrixMulTri@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatMulTri_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE
