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

!>  \brief This file contains the implementations of the tests of module [pm_matrixMulAdd](@ref pm_matrixMulAdd).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_matrixMulAdd) routines

    use pm_io, only: display_type
    use pm_mathCompare, only: isClose
    use pm_matrixCopy, only: rdpack
    use pm_matrixCopy, only: getMatCopy
    use pm_matrixCopy, only: setMatCopy
    use pm_arrayChoice, only: getChoice
    use pm_arrayResize, only: setResized
    use pm_arrayRebill, only: setRebilled
    use pm_distUnif, only: getUnifRand, setUnifRand

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatMulAdd_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setMatMulAdd_IK5
        use pm_kind, only: IKC => IK5
#include "test_pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setMatMulAdd_IK4
        use pm_kind, only: IKC => IK4
#include "test_pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setMatMulAdd_IK3
        use pm_kind, only: IKC => IK3
#include "test_pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setMatMulAdd_IK2
        use pm_kind, only: IKC => IK2
#include "test_pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setMatMulAdd_IK1
        use pm_kind, only: IKC => IK1
#include "test_pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setMatMulAdd_CK5
        use pm_kind, only: CKC => CK5
#include "test_pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setMatMulAdd_CK4
        use pm_kind, only: CKC => CK4
#include "test_pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setMatMulAdd_CK3
        use pm_kind, only: CKC => CK3
#include "test_pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setMatMulAdd_CK2
        use pm_kind, only: CKC => CK2
#include "test_pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setMatMulAdd_CK1
        use pm_kind, only: CKC => CK1
#include "test_pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setMatMulAdd_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setMatMulAdd_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setMatMulAdd_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setMatMulAdd_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setMatMulAdd_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_matrixMulAdd@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatMulAdd_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE
