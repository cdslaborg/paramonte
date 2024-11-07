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
!>  This module contains tests of the module [pm_arrayUnique](@ref pm_arrayUnique).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_arrayUnique

    use pm_arrayUnique
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_isUnique_D0_SK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_isUnique_D0_SK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_isUnique_D0_SK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_isUnique_D0_SK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_isUnique_D0_SK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_isUnique_D1_SK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_isUnique_D1_SK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_isUnique_D1_SK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_isUnique_D1_SK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_isUnique_D1_SK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_isUnique_D1_IK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_isUnique_D1_IK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_isUnique_D1_IK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_isUnique_D1_IK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_isUnique_D1_IK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_isUnique_D1_LK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_isUnique_D1_LK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_isUnique_D1_LK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_isUnique_D1_LK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_isUnique_D1_LK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_isUnique_D1_CK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_isUnique_D1_CK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_isUnique_D1_CK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_isUnique_D1_CK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_isUnique_D1_CK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_isUnique_D1_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_isUnique_D1_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_isUnique_D1_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_isUnique_D1_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_isUnique_D1_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_isUniqueAll_D0_SK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_isUniqueAll_D0_SK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_isUniqueAll_D0_SK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_isUniqueAll_D0_SK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_isUniqueAll_D0_SK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_isUniqueAll_D1_SK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_isUniqueAll_D1_SK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_isUniqueAll_D1_SK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_isUniqueAll_D1_SK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_isUniqueAll_D1_SK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_isUniqueAll_D1_IK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_isUniqueAll_D1_IK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_isUniqueAll_D1_IK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_isUniqueAll_D1_IK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_isUniqueAll_D1_IK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_isUniqueAll_D1_LK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_isUniqueAll_D1_LK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_isUniqueAll_D1_LK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_isUniqueAll_D1_LK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_isUniqueAll_D1_LK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_isUniqueAll_D1_CK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_isUniqueAll_D1_CK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_isUniqueAll_D1_CK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_isUniqueAll_D1_CK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_isUniqueAll_D1_CK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_isUniqueAll_D1_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_isUniqueAll_D1_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_isUniqueAll_D1_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_isUniqueAll_D1_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_isUniqueAll_D1_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_isUniqueAny_D0_SK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_isUniqueAny_D0_SK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_isUniqueAny_D0_SK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_isUniqueAny_D0_SK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_isUniqueAny_D0_SK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_isUniqueAny_D1_SK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_isUniqueAny_D1_SK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_isUniqueAny_D1_SK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_isUniqueAny_D1_SK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_isUniqueAny_D1_SK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_isUniqueAny_D1_IK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_isUniqueAny_D1_IK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_isUniqueAny_D1_IK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_isUniqueAny_D1_IK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_isUniqueAny_D1_IK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_isUniqueAny_D1_LK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_isUniqueAny_D1_LK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_isUniqueAny_D1_LK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_isUniqueAny_D1_LK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_isUniqueAny_D1_LK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_isUniqueAny_D1_CK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_isUniqueAny_D1_CK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_isUniqueAny_D1_CK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_isUniqueAny_D1_CK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_isUniqueAny_D1_CK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_isUniqueAny_D1_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_isUniqueAny_D1_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_isUniqueAny_D1_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_isUniqueAny_D1_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_isUniqueAny_D1_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_getUnique_D0_SK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_getUnique_D0_SK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_getUnique_D0_SK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_getUnique_D0_SK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_getUnique_D0_SK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_getUnique_D1_SK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_getUnique_D1_SK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_getUnique_D1_SK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_getUnique_D1_SK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_getUnique_D1_SK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_getUnique_D1_IK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_getUnique_D1_IK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_getUnique_D1_IK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_getUnique_D1_IK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_getUnique_D1_IK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_getUnique_D1_LK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_getUnique_D1_LK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_getUnique_D1_LK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_getUnique_D1_LK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_getUnique_D1_LK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_getUnique_D1_CK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_getUnique_D1_CK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_getUnique_D1_CK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_getUnique_D1_CK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_getUnique_D1_CK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getUnique_D1_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getUnique_D1_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getUnique_D1_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getUnique_D1_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getUnique_D1_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setUnique_D0_SK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setUnique_D0_SK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setUnique_D0_SK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setUnique_D0_SK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setUnique_D0_SK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setUnique_D1_SK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setUnique_D1_SK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setUnique_D1_SK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setUnique_D1_SK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setUnique_D1_SK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setUnique_D1_IK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setUnique_D1_IK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setUnique_D1_IK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setUnique_D1_IK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setUnique_D1_IK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_setUnique_D1_LK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setUnique_D1_LK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setUnique_D1_LK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setUnique_D1_LK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setUnique_D1_LK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setUnique_D1_CK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setUnique_D1_CK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setUnique_D1_CK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setUnique_D1_CK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setUnique_D1_CK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setUnique_D1_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setUnique_D1_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setUnique_D1_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setUnique_D1_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setUnique_D1_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none

        test = test_type(MODULE_NAME)

        call test%run(test_getUnique_D1_1, SK_"test_getUnique_D1_1")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_isUnique_D0_SK5_1, SK_"test_isUnique_D0_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_isUnique_D0_SK4_1, SK_"test_isUnique_D0_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_isUnique_D0_SK3_1, SK_"test_isUnique_D0_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_isUnique_D0_SK2_1, SK_"test_isUnique_D0_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_isUnique_D0_SK1_1, SK_"test_isUnique_D0_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_isUnique_D1_SK5_1, SK_"test_isUnique_D1_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_isUnique_D1_SK4_1, SK_"test_isUnique_D1_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_isUnique_D1_SK3_1, SK_"test_isUnique_D1_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_isUnique_D1_SK2_1, SK_"test_isUnique_D1_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_isUnique_D1_SK1_1, SK_"test_isUnique_D1_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_isUnique_D1_IK5_1, SK_"test_isUnique_D1_IK5_1")
#endif
#if     IK4_ENABLED
        call test%run(test_isUnique_D1_IK4_1, SK_"test_isUnique_D1_IK4_1")
#endif
#if     IK3_ENABLED
        call test%run(test_isUnique_D1_IK3_1, SK_"test_isUnique_D1_IK3_1")
#endif
#if     IK2_ENABLED
        call test%run(test_isUnique_D1_IK2_1, SK_"test_isUnique_D1_IK2_1")
#endif
#if     IK1_ENABLED
        call test%run(test_isUnique_D1_IK1_1, SK_"test_isUnique_D1_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_isUnique_D1_LK5_1, SK_"test_isUnique_D1_LK5_1")
#endif
#if     LK4_ENABLED
        call test%run(test_isUnique_D1_LK4_1, SK_"test_isUnique_D1_LK4_1")
#endif
#if     LK3_ENABLED
        call test%run(test_isUnique_D1_LK3_1, SK_"test_isUnique_D1_LK3_1")
#endif
#if     LK2_ENABLED
        call test%run(test_isUnique_D1_LK2_1, SK_"test_isUnique_D1_LK2_1")
#endif
#if     LK1_ENABLED
        call test%run(test_isUnique_D1_LK1_1, SK_"test_isUnique_D1_LK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_isUnique_D1_CK5_1, SK_"test_isUnique_D1_CK5_1")
#endif
#if     CK4_ENABLED
        call test%run(test_isUnique_D1_CK4_1, SK_"test_isUnique_D1_CK4_1")
#endif
#if     CK3_ENABLED
        call test%run(test_isUnique_D1_CK3_1, SK_"test_isUnique_D1_CK3_1")
#endif
#if     CK2_ENABLED
        call test%run(test_isUnique_D1_CK2_1, SK_"test_isUnique_D1_CK2_1")
#endif
#if     CK1_ENABLED
        call test%run(test_isUnique_D1_CK1_1, SK_"test_isUnique_D1_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_isUnique_D1_RK5_1, SK_"test_isUnique_D1_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_isUnique_D1_RK4_1, SK_"test_isUnique_D1_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_isUnique_D1_RK3_1, SK_"test_isUnique_D1_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_isUnique_D1_RK2_1, SK_"test_isUnique_D1_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_isUnique_D1_RK1_1, SK_"test_isUnique_D1_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_isUniqueAll_D0_SK5_1, SK_"test_isUniqueAll_D0_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_isUniqueAll_D0_SK4_1, SK_"test_isUniqueAll_D0_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_isUniqueAll_D0_SK3_1, SK_"test_isUniqueAll_D0_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_isUniqueAll_D0_SK2_1, SK_"test_isUniqueAll_D0_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_isUniqueAll_D0_SK1_1, SK_"test_isUniqueAll_D0_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_isUniqueAll_D1_SK5_1, SK_"test_isUniqueAll_D1_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_isUniqueAll_D1_SK4_1, SK_"test_isUniqueAll_D1_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_isUniqueAll_D1_SK3_1, SK_"test_isUniqueAll_D1_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_isUniqueAll_D1_SK2_1, SK_"test_isUniqueAll_D1_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_isUniqueAll_D1_SK1_1, SK_"test_isUniqueAll_D1_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_isUniqueAll_D1_IK5_1, SK_"test_isUniqueAll_D1_IK5_1")
#endif
#if     IK4_ENABLED
        call test%run(test_isUniqueAll_D1_IK4_1, SK_"test_isUniqueAll_D1_IK4_1")
#endif
#if     IK3_ENABLED
        call test%run(test_isUniqueAll_D1_IK3_1, SK_"test_isUniqueAll_D1_IK3_1")
#endif
#if     IK2_ENABLED
        call test%run(test_isUniqueAll_D1_IK2_1, SK_"test_isUniqueAll_D1_IK2_1")
#endif
#if     IK1_ENABLED
        call test%run(test_isUniqueAll_D1_IK1_1, SK_"test_isUniqueAll_D1_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_isUniqueAll_D1_LK5_1, SK_"test_isUniqueAll_D1_LK5_1")
#endif
#if     LK4_ENABLED
        call test%run(test_isUniqueAll_D1_LK4_1, SK_"test_isUniqueAll_D1_LK4_1")
#endif
#if     LK3_ENABLED
        call test%run(test_isUniqueAll_D1_LK3_1, SK_"test_isUniqueAll_D1_LK3_1")
#endif
#if     LK2_ENABLED
        call test%run(test_isUniqueAll_D1_LK2_1, SK_"test_isUniqueAll_D1_LK2_1")
#endif
#if     LK1_ENABLED
        call test%run(test_isUniqueAll_D1_LK1_1, SK_"test_isUniqueAll_D1_LK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_isUniqueAll_D1_CK5_1, SK_"test_isUniqueAll_D1_CK5_1")
#endif
#if     CK4_ENABLED
        call test%run(test_isUniqueAll_D1_CK4_1, SK_"test_isUniqueAll_D1_CK4_1")
#endif
#if     CK3_ENABLED
        call test%run(test_isUniqueAll_D1_CK3_1, SK_"test_isUniqueAll_D1_CK3_1")
#endif
#if     CK2_ENABLED
        call test%run(test_isUniqueAll_D1_CK2_1, SK_"test_isUniqueAll_D1_CK2_1")
#endif
#if     CK1_ENABLED
        call test%run(test_isUniqueAll_D1_CK1_1, SK_"test_isUniqueAll_D1_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_isUniqueAll_D1_RK5_1, SK_"test_isUniqueAll_D1_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_isUniqueAll_D1_RK4_1, SK_"test_isUniqueAll_D1_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_isUniqueAll_D1_RK3_1, SK_"test_isUniqueAll_D1_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_isUniqueAll_D1_RK2_1, SK_"test_isUniqueAll_D1_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_isUniqueAll_D1_RK1_1, SK_"test_isUniqueAll_D1_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_isUniqueAny_D0_SK5_1, SK_"test_isUniqueAny_D0_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_isUniqueAny_D0_SK4_1, SK_"test_isUniqueAny_D0_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_isUniqueAny_D0_SK3_1, SK_"test_isUniqueAny_D0_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_isUniqueAny_D0_SK2_1, SK_"test_isUniqueAny_D0_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_isUniqueAny_D0_SK1_1, SK_"test_isUniqueAny_D0_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_isUniqueAny_D1_SK5_1, SK_"test_isUniqueAny_D1_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_isUniqueAny_D1_SK4_1, SK_"test_isUniqueAny_D1_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_isUniqueAny_D1_SK3_1, SK_"test_isUniqueAny_D1_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_isUniqueAny_D1_SK2_1, SK_"test_isUniqueAny_D1_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_isUniqueAny_D1_SK1_1, SK_"test_isUniqueAny_D1_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_isUniqueAny_D1_IK5_1, SK_"test_isUniqueAny_D1_IK5_1")
#endif
#if     IK4_ENABLED
        call test%run(test_isUniqueAny_D1_IK4_1, SK_"test_isUniqueAny_D1_IK4_1")
#endif
#if     IK3_ENABLED
        call test%run(test_isUniqueAny_D1_IK3_1, SK_"test_isUniqueAny_D1_IK3_1")
#endif
#if     IK2_ENABLED
        call test%run(test_isUniqueAny_D1_IK2_1, SK_"test_isUniqueAny_D1_IK2_1")
#endif
#if     IK1_ENABLED
        call test%run(test_isUniqueAny_D1_IK1_1, SK_"test_isUniqueAny_D1_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_isUniqueAny_D1_LK5_1, SK_"test_isUniqueAny_D1_LK5_1")
#endif
#if     LK4_ENABLED
        call test%run(test_isUniqueAny_D1_LK4_1, SK_"test_isUniqueAny_D1_LK4_1")
#endif
#if     LK3_ENABLED
        call test%run(test_isUniqueAny_D1_LK3_1, SK_"test_isUniqueAny_D1_LK3_1")
#endif
#if     LK2_ENABLED
        call test%run(test_isUniqueAny_D1_LK2_1, SK_"test_isUniqueAny_D1_LK2_1")
#endif
#if     LK1_ENABLED
        call test%run(test_isUniqueAny_D1_LK1_1, SK_"test_isUniqueAny_D1_LK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_isUniqueAny_D1_CK5_1, SK_"test_isUniqueAny_D1_CK5_1")
#endif
#if     CK4_ENABLED
        call test%run(test_isUniqueAny_D1_CK4_1, SK_"test_isUniqueAny_D1_CK4_1")
#endif
#if     CK3_ENABLED
        call test%run(test_isUniqueAny_D1_CK3_1, SK_"test_isUniqueAny_D1_CK3_1")
#endif
#if     CK2_ENABLED
        call test%run(test_isUniqueAny_D1_CK2_1, SK_"test_isUniqueAny_D1_CK2_1")
#endif
#if     CK1_ENABLED
        call test%run(test_isUniqueAny_D1_CK1_1, SK_"test_isUniqueAny_D1_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_isUniqueAny_D1_RK5_1, SK_"test_isUniqueAny_D1_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_isUniqueAny_D1_RK4_1, SK_"test_isUniqueAny_D1_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_isUniqueAny_D1_RK3_1, SK_"test_isUniqueAny_D1_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_isUniqueAny_D1_RK2_1, SK_"test_isUniqueAny_D1_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_isUniqueAny_D1_RK1_1, SK_"test_isUniqueAny_D1_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_getUnique_D0_SK5_1, SK_"test_getUnique_D0_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_getUnique_D0_SK4_1, SK_"test_getUnique_D0_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_getUnique_D0_SK3_1, SK_"test_getUnique_D0_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_getUnique_D0_SK2_1, SK_"test_getUnique_D0_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_getUnique_D0_SK1_1, SK_"test_getUnique_D0_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_getUnique_D1_SK5_1, SK_"test_getUnique_D1_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_getUnique_D1_SK4_1, SK_"test_getUnique_D1_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_getUnique_D1_SK3_1, SK_"test_getUnique_D1_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_getUnique_D1_SK2_1, SK_"test_getUnique_D1_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_getUnique_D1_SK1_1, SK_"test_getUnique_D1_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_getUnique_D1_IK5_1, SK_"test_getUnique_D1_IK5_1")
#endif
#if     IK4_ENABLED
        call test%run(test_getUnique_D1_IK4_1, SK_"test_getUnique_D1_IK4_1")
#endif
#if     IK3_ENABLED
        call test%run(test_getUnique_D1_IK3_1, SK_"test_getUnique_D1_IK3_1")
#endif
#if     IK2_ENABLED
        call test%run(test_getUnique_D1_IK2_1, SK_"test_getUnique_D1_IK2_1")
#endif
#if     IK1_ENABLED
        call test%run(test_getUnique_D1_IK1_1, SK_"test_getUnique_D1_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_getUnique_D1_LK5_1, SK_"test_getUnique_D1_LK5_1")
#endif
#if     LK4_ENABLED
        call test%run(test_getUnique_D1_LK4_1, SK_"test_getUnique_D1_LK4_1")
#endif
#if     LK3_ENABLED
        call test%run(test_getUnique_D1_LK3_1, SK_"test_getUnique_D1_LK3_1")
#endif
#if     LK2_ENABLED
        call test%run(test_getUnique_D1_LK2_1, SK_"test_getUnique_D1_LK2_1")
#endif
#if     LK1_ENABLED
        call test%run(test_getUnique_D1_LK1_1, SK_"test_getUnique_D1_LK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_getUnique_D1_CK5_1, SK_"test_getUnique_D1_CK5_1")
#endif
#if     CK4_ENABLED
        call test%run(test_getUnique_D1_CK4_1, SK_"test_getUnique_D1_CK4_1")
#endif
#if     CK3_ENABLED
        call test%run(test_getUnique_D1_CK3_1, SK_"test_getUnique_D1_CK3_1")
#endif
#if     CK2_ENABLED
        call test%run(test_getUnique_D1_CK2_1, SK_"test_getUnique_D1_CK2_1")
#endif
#if     CK1_ENABLED
        call test%run(test_getUnique_D1_CK1_1, SK_"test_getUnique_D1_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getUnique_D1_RK5_1, SK_"test_getUnique_D1_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_getUnique_D1_RK4_1, SK_"test_getUnique_D1_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_getUnique_D1_RK3_1, SK_"test_getUnique_D1_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_getUnique_D1_RK2_1, SK_"test_getUnique_D1_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_getUnique_D1_RK1_1, SK_"test_getUnique_D1_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setUnique_D0_SK5_1, SK_"test_setUnique_D0_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_setUnique_D0_SK4_1, SK_"test_setUnique_D0_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_setUnique_D0_SK3_1, SK_"test_setUnique_D0_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_setUnique_D0_SK2_1, SK_"test_setUnique_D0_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_setUnique_D0_SK1_1, SK_"test_setUnique_D0_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setUnique_D1_SK5_1, SK_"test_setUnique_D1_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_setUnique_D1_SK4_1, SK_"test_setUnique_D1_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_setUnique_D1_SK3_1, SK_"test_setUnique_D1_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_setUnique_D1_SK2_1, SK_"test_setUnique_D1_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_setUnique_D1_SK1_1, SK_"test_setUnique_D1_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setUnique_D1_IK5_1, SK_"test_setUnique_D1_IK5_1")
#endif
#if     IK4_ENABLED
        call test%run(test_setUnique_D1_IK4_1, SK_"test_setUnique_D1_IK4_1")
#endif
#if     IK3_ENABLED
        call test%run(test_setUnique_D1_IK3_1, SK_"test_setUnique_D1_IK3_1")
#endif
#if     IK2_ENABLED
        call test%run(test_setUnique_D1_IK2_1, SK_"test_setUnique_D1_IK2_1")
#endif
#if     IK1_ENABLED
        call test%run(test_setUnique_D1_IK1_1, SK_"test_setUnique_D1_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setUnique_D1_LK5_1, SK_"test_setUnique_D1_LK5_1")
#endif
#if     LK4_ENABLED
        call test%run(test_setUnique_D1_LK4_1, SK_"test_setUnique_D1_LK4_1")
#endif
#if     LK3_ENABLED
        call test%run(test_setUnique_D1_LK3_1, SK_"test_setUnique_D1_LK3_1")
#endif
#if     LK2_ENABLED
        call test%run(test_setUnique_D1_LK2_1, SK_"test_setUnique_D1_LK2_1")
#endif
#if     LK1_ENABLED
        call test%run(test_setUnique_D1_LK1_1, SK_"test_setUnique_D1_LK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setUnique_D1_CK5_1, SK_"test_setUnique_D1_CK5_1")
#endif
#if     CK4_ENABLED
        call test%run(test_setUnique_D1_CK4_1, SK_"test_setUnique_D1_CK4_1")
#endif
#if     CK3_ENABLED
        call test%run(test_setUnique_D1_CK3_1, SK_"test_setUnique_D1_CK3_1")
#endif
#if     CK2_ENABLED
        call test%run(test_setUnique_D1_CK2_1, SK_"test_setUnique_D1_CK2_1")
#endif
#if     CK1_ENABLED
        call test%run(test_setUnique_D1_CK1_1, SK_"test_setUnique_D1_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setUnique_D1_RK5_1, SK_"test_setUnique_D1_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setUnique_D1_RK4_1, SK_"test_setUnique_D1_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setUnique_D1_RK3_1, SK_"test_setUnique_D1_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setUnique_D1_RK2_1, SK_"test_setUnique_D1_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setUnique_D1_RK1_1, SK_"test_setUnique_D1_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getUnique_D1_1() result(assertion)

        use pm_kind, only: RK, IK

        implicit none
        logical(LK)                 :: assertion
        integer(IK) , parameter     :: VECTOR(*) = int([1,2,1,3,5,5,2],IK)
        integer(IK) , parameter     :: UNIQUE_VALUE(*) = int([1,2,3,5],IK)
        integer(IK) , allocatable   :: Unique(:)

        Unique = getUnique(VECTOR)
        assertion = all(Unique == UNIQUE_VALUE)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "VECTOR", VECTOR
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "UNIQUE_VALUE   ", UNIQUE_VALUE
            write(test%disp%unit,"(*(g0,:,', '))") "Unique         ", Unique
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getUnique_D1_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !function test_findUnique_1() result(assertion)
    !
    !    use pm_kind, only: RK, IK
    !    use pm_container, only: IV => cvi_pdt
    !
    !    implicit none
    !    logical(LK)                 :: assertion
    !    logical(LK)                 :: assertionCurrent
    !    integer(IK) , parameter     :: VECTOR(*) = [1,2,1,3,5,5,2]
    !    integer(IK) , parameter     :: LEN_VECTOR = size(VECTOR)
    !    integer(IK) , parameter     :: UNIQUE_VALUE(*) = [1,2,3,5]
    !    integer(IK) , parameter     :: UNIQUE_COUNT(*) = [2,2,1,2]
    !    integer(IK) , allocatable   :: ZeroLenVector(:)
    !    integer(IK) , allocatable   :: UniqueValue(:)
    !    integer(IK) , allocatable   :: UniqueCount(:)
    !    type(IV)    , allocatable   :: UniqueIndex(:)
    !    type(IV)    , allocatable   :: UNIQUE_INDEX(:)
    !    integer(IK)                 :: lenUnique, i
    !    type(err_type)              :: Err
    !
    !    call findUnique ( Vector = VECTOR &
    !                    , lenUnique = lenUnique &
    !                    , UniqueValue = UniqueValue &
    !                    , UniqueCount = UniqueCount &
    !                    )
    !
    !    assertion = all(UniqueValue(1:lenUnique)==UNIQUE_VALUE) .and. all(UniqueCount(1:lenUnique)==UNIQUE_COUNT)
    !
    !    if (test%traceable .and. .not. assertion) then
    !    ! LCOV_EXCL_START
    !        write(test%disp%unit,"(*(g0,:,', '))")
    !        write(test%disp%unit,"(*(g0,:,', '))") "VECTOR", VECTOR
    !        write(test%disp%unit,"(*(g0,:,', '))")
    !        write(test%disp%unit,"(*(g0,:,', '))") "UNIQUE_VALUE", UNIQUE_VALUE
    !        write(test%disp%unit,"(*(g0,:,', '))") "UniqueValue ", UniqueValue(1:lenUnique)
    !        write(test%disp%unit,"(*(g0,:,', '))")
    !        write(test%disp%unit,"(*(g0,:,', '))") "UNIQUE_COUNT", UNIQUE_COUNT
    !        write(test%disp%unit,"(*(g0,:,', '))") "UniqueCount ", UniqueCount(1:lenUnique)
    !        write(test%disp%unit,"(*(g0,:,', '))")
    !        write(test%disp%unit,"(*(g0,:,', '))") "lenUnique", lenUnique
    !        write(test%disp%unit,"(*(g0,:,', '))")
    !    end if
    !    ! LCOV_EXCL_STOP
    !
    !    call test%assert(assertion)
    !
    !    ! test UniqueIndex
    !
    !    allocate(UNIQUE_INDEX(size(UNIQUE_COUNT)))
    !    UNIQUE_INDEX(1)%val = [5,6]
    !    UNIQUE_INDEX(2)%val = [2,7]
    !    UNIQUE_INDEX(3)%val = [1,3]
    !    UNIQUE_INDEX(4)%val = [4]
    !
    !    call findUnique ( Vector = VECTOR &
    !                    , lenUnique = lenUnique &
    !                    , UniqueValue = UniqueValue &
    !                    , UniqueCount = UniqueCount &
    !                    , UniqueIndex = UniqueIndex &
    !                    , sorting = -1_IK &
    !                    )
    !
    !    assertion = assertion .and. .not. err%occurred
    !    call test%assert(assertion)
    !
    !    do i = 1, lenUnique
    !
    !        assertionCurrent = all(UniqueIndex(i)%val == UNIQUE_INDEX(i)%val)
    !        assertion = assertion .and. assertionCurrent
    !
    !        if (test%traceable .and. .not. assertionCurrent) then
    !        ! LCOV_EXCL_START
    !            write(test%disp%unit,"(*(g0,:,', '))")
    !            write(test%disp%unit,"(*(g0,:,', '))") "UNIQUE_COUNT", UNIQUE_COUNT
    !            write(test%disp%unit,"(*(g0,:,', '))") "UniqueCount ", UniqueCount
    !            write(test%disp%unit,"(*(g0,:,', '))")
    !        end if
    !        ! LCOV_EXCL_STOP
    !
    !        if (i>1_IK) assertionCurrent = assertionCurrent .and. UniqueCount(i) <= UniqueCount(i-1)
    !        assertion = assertion .and. assertionCurrent
    !
    !        if (test%traceable .and. .not. assertionCurrent) then
    !        ! LCOV_EXCL_START
    !            write(test%disp%unit,"(*(g0,:,', '))")
    !            write(test%disp%unit,"(*(g0,:,', '))") "VECTOR", VECTOR
    !            write(test%disp%unit,"(*(g0,:,', '))")
    !            write(test%disp%unit,"(*(g0,:,', '))") "UNIQUE_VALUE", UNIQUE_VALUE
    !            write(test%disp%unit,"(*(g0,:,', '))") "UniqueValue ", UniqueValue
    !            write(test%disp%unit,"(*(g0,:,', '))")
    !            write(test%disp%unit,"(*(g0,:,', '))") "UNIQUE_COUNT", UNIQUE_COUNT
    !            write(test%disp%unit,"(*(g0,:,', '))") "UniqueCount ", UniqueCount
    !            write(test%disp%unit,"(*(g0,:,', '))")
    !            write(test%disp%unit,"(*(g0,:,', '))") "UNIQUE_INDEX(i)%val  ", UNIQUE_INDEX(i)%val
    !            write(test%disp%unit,"(*(g0,:,', '))") "UniqueIndex(i)%val   ", UniqueIndex(i)%val
    !            write(test%disp%unit,"(*(g0,:,', '))")
    !            write(test%disp%unit,"(*(g0,:,', '))") "lenUnique", lenUnique
    !            write(test%disp%unit,"(*(g0,:,', '))")
    !        end if
    !        ! LCOV_EXCL_STOP
    !
    !        if (.not. assertion) exit ! LCOV_EXCL_LINE
    !
    !    end do
    !    call test%assert(assertion)
    !
    !    ! test with empty input vector
    !
    !    allocate(ZeroLenVector(0))
    !    call findUnique ( Vector = ZeroLenVector & ! LCOV_EXCL_LINE
    !                    , UniqueValue = UniqueValue & ! LCOV_EXCL_LINE
    !                    , UniqueCount = UniqueCount & ! LCOV_EXCL_LINE
    !                    , lenUnique = lenUnique & ! LCOV_EXCL_LINE
    !                    )
    !
    !    if (test%traceable .and. .not. assertion) then
    !    ! LCOV_EXCL_START
    !        write(test%disp%unit,"(*(g0,:,', '))")
    !        write(test%disp%unit,"(*(g0,:,', '))") "VECTOR", VECTOR
    !        write(test%disp%unit,"(*(g0,:,', '))")
    !        write(test%disp%unit,"(*(g0,:,', '))") "UNIQUE_VALUE", UNIQUE_VALUE
    !        write(test%disp%unit,"(*(g0,:,', '))") "UniqueValue ", UniqueValue(1:lenUnique)
    !        write(test%disp%unit,"(*(g0,:,', '))")
    !        write(test%disp%unit,"(*(g0,:,', '))") "UNIQUE_COUNT", UNIQUE_COUNT
    !        write(test%disp%unit,"(*(g0,:,', '))") "UniqueCount ", UniqueCount(1:lenUnique)
    !        write(test%disp%unit,"(*(g0,:,', '))")
    !        write(test%disp%unit,"(*(g0,:,', '))") "lenUnique", lenUnique
    !        write(test%disp%unit,"(*(g0,:,', '))")
    !        write(test%disp%unit,"(*(g0,:,', '))")
    !        write(test%disp%unit,"(*(g0,:,', '))")
    !        write(test%disp%unit,"(*(g0,:,', '))") "VECTOR", ZeroLenVector
    !        write(test%disp%unit,"(*(g0,:,', '))")
    !        write(test%disp%unit,"(*(g0,:,', '))") "UniqueValue ", UniqueValue(1:lenUnique)
    !        write(test%disp%unit,"(*(g0,:,', '))")
    !        write(test%disp%unit,"(*(g0,:,', '))") "UniqueCount ", UniqueCount(1:lenUnique)
    !        write(test%disp%unit,"(*(g0,:,', '))")
    !        write(test%disp%unit,"(*(g0,:,', '))") "lenUnique", lenUnique
    !        write(test%disp%unit,"(*(g0,:,', '))")
    !    end if
    !    ! LCOV_EXCL_STOP
    !
    !end function test_findUnique_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_arrayUnique ! LCOV_EXCL_LINE