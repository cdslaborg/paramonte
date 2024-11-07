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

!>  \brief This file contains the implementations of the tests of module [pm_distanceEuclid](@ref pm_distanceEuclid).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, March 22, 2012, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_distanceEuclid) routines

    use pm_container, only: csp_type, css_type
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: getUnifRand
    use pm_arrayFill, only: getFilled
    use pm_arrayRange, only: getRange

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisEuclid_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getDisEuclid_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getDisEuclid_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getDisEuclid_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getDisEuclid_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getDisEuclid_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisEuclid_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setDisEuclid_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setDisEuclid_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setDisEuclid_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setDisEuclid_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setDisEuclid_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setDisEuclid_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setDisEuclid_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisMatEuclid_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getDisMatEuclid_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getDisMatEuclid_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getDisMatEuclid_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getDisMatEuclid_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getDisMatEuclid_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisMatEuclid_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setDisMatEuclid_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setDisMatEuclid_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setDisMatEuclid_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setDisMatEuclid_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setDisMatEuclid_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setDisMatEuclid_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setDisMatEuclid_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE