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
!>  This module contains tests of the module [pm_arrayReplace](@ref pm_arrayReplace).
!>
!>  \todo
!>  \phigh The following tests should be extended to include allocatable arrays of arbitrary bounds.
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_arrayReplace

    use pm_arrayReplace
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_getReplaced_D0_D0_D0_SK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D0_D0_D0_SK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_getReplaced_D0_D0_D0_SK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D0_D0_D0_SK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_getReplaced_D0_D0_D0_SK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D0_D0_D0_SK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_getReplaced_D0_D0_D0_SK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D0_D0_D0_SK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_getReplaced_D0_D0_D0_SK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D0_D0_D0_SK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_getReplaced_D1_D1_D1_SK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_SK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_getReplaced_D1_D1_D1_SK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_SK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_getReplaced_D1_D1_D1_SK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_SK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_getReplaced_D1_D1_D1_SK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_SK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_getReplaced_D1_D1_D1_SK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_SK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_getReplaced_D1_D1_D1_IK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_IK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_getReplaced_D1_D1_D1_IK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_IK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_getReplaced_D1_D1_D1_IK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_IK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_getReplaced_D1_D1_D1_IK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_IK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_getReplaced_D1_D1_D1_IK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_IK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_getReplaced_D1_D1_D1_LK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_LK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_getReplaced_D1_D1_D1_LK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_LK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_getReplaced_D1_D1_D1_LK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_LK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_getReplaced_D1_D1_D1_LK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_LK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_getReplaced_D1_D1_D1_LK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_LK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_getReplaced_D1_D1_D1_CK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_CK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_getReplaced_D1_D1_D1_CK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_CK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_getReplaced_D1_D1_D1_CK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_CK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_getReplaced_D1_D1_D1_CK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_CK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_getReplaced_D1_D1_D1_CK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_CK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getReplaced_D1_D1_D1_RK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_RK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getReplaced_D1_D1_D1_RK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_RK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getReplaced_D1_D1_D1_RK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_RK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getReplaced_D1_D1_D1_RK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_RK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getReplaced_D1_D1_D1_RK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D1_RK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_getReplaced_D1_D0_D0_SK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_SK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_getReplaced_D1_D0_D0_SK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_SK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_getReplaced_D1_D0_D0_SK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_SK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_getReplaced_D1_D0_D0_SK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_SK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_getReplaced_D1_D0_D0_SK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_SK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_getReplaced_D1_D0_D0_IK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_IK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_getReplaced_D1_D0_D0_IK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_IK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_getReplaced_D1_D0_D0_IK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_IK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_getReplaced_D1_D0_D0_IK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_IK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_getReplaced_D1_D0_D0_IK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_IK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_getReplaced_D1_D0_D0_LK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_LK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_getReplaced_D1_D0_D0_LK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_LK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_getReplaced_D1_D0_D0_LK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_LK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_getReplaced_D1_D0_D0_LK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_LK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_getReplaced_D1_D0_D0_LK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_LK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_getReplaced_D1_D0_D0_CK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_CK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_getReplaced_D1_D0_D0_CK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_CK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_getReplaced_D1_D0_D0_CK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_CK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_getReplaced_D1_D0_D0_CK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_CK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_getReplaced_D1_D0_D0_CK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_CK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getReplaced_D1_D0_D0_RK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_RK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getReplaced_D1_D0_D0_RK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_RK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getReplaced_D1_D0_D0_RK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_RK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getReplaced_D1_D0_D0_RK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_RK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getReplaced_D1_D0_D0_RK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D0_RK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_getReplaced_D1_D0_D1_SK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_SK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_getReplaced_D1_D0_D1_SK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_SK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_getReplaced_D1_D0_D1_SK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_SK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_getReplaced_D1_D0_D1_SK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_SK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_getReplaced_D1_D0_D1_SK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_SK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_getReplaced_D1_D0_D1_IK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_IK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_getReplaced_D1_D0_D1_IK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_IK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_getReplaced_D1_D0_D1_IK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_IK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_getReplaced_D1_D0_D1_IK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_IK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_getReplaced_D1_D0_D1_IK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_IK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_getReplaced_D1_D0_D1_LK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_LK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_getReplaced_D1_D0_D1_LK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_LK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_getReplaced_D1_D0_D1_LK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_LK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_getReplaced_D1_D0_D1_LK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_LK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_getReplaced_D1_D0_D1_LK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_LK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_getReplaced_D1_D0_D1_CK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_CK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_getReplaced_D1_D0_D1_CK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_CK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_getReplaced_D1_D0_D1_CK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_CK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_getReplaced_D1_D0_D1_CK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_CK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_getReplaced_D1_D0_D1_CK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_CK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getReplaced_D1_D0_D1_RK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_RK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getReplaced_D1_D0_D1_RK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_RK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getReplaced_D1_D0_D1_RK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_RK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getReplaced_D1_D0_D1_RK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_RK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getReplaced_D1_D0_D1_RK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D0_D1_RK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_getReplaced_D1_D1_D0_SK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_SK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_getReplaced_D1_D1_D0_SK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_SK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_getReplaced_D1_D1_D0_SK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_SK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_getReplaced_D1_D1_D0_SK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_SK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_getReplaced_D1_D1_D0_SK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_SK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_getReplaced_D1_D1_D0_IK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_IK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_getReplaced_D1_D1_D0_IK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_IK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_getReplaced_D1_D1_D0_IK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_IK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_getReplaced_D1_D1_D0_IK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_IK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_getReplaced_D1_D1_D0_IK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_IK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_getReplaced_D1_D1_D0_LK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_LK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_getReplaced_D1_D1_D0_LK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_LK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_getReplaced_D1_D1_D0_LK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_LK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_getReplaced_D1_D1_D0_LK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_LK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_getReplaced_D1_D1_D0_LK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_LK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_getReplaced_D1_D1_D0_CK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_CK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_getReplaced_D1_D1_D0_CK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_CK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_getReplaced_D1_D1_D0_CK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_CK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_getReplaced_D1_D1_D0_CK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_CK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_getReplaced_D1_D1_D0_CK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_CK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getReplaced_D1_D1_D0_RK5_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_RK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getReplaced_D1_D1_D0_RK4_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_RK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getReplaced_D1_D1_D0_RK3_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_RK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getReplaced_D1_D1_D0_RK2_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_RK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getReplaced_D1_D1_D0_RK1_1 () result(assertion); logical(LK) :: assertion; end function
        module function test_setReplaced_D1_D1_D0_RK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        !call test%run(test_getReplaced_D0_D0_D0_SK5_1, SK_"test_getReplaced_D0_D0_D0_SK5_1")
        call test%run(test_setReplaced_D0_D0_D0_SK5_1, SK_"test_setReplaced_D0_D0_D0_SK5_1")
#endif
#if     SK4_ENABLED
        !call test%run(test_getReplaced_D0_D0_D0_SK4_1, SK_"test_getReplaced_D0_D0_D0_SK4_1")
        call test%run(test_setReplaced_D0_D0_D0_SK4_1, SK_"test_setReplaced_D0_D0_D0_SK4_1")
#endif
#if     SK3_ENABLED
        !call test%run(test_getReplaced_D0_D0_D0_SK3_1, SK_"test_getReplaced_D0_D0_D0_SK3_1")
        call test%run(test_setReplaced_D0_D0_D0_SK3_1, SK_"test_setReplaced_D0_D0_D0_SK3_1")
#endif
#if     SK2_ENABLED
        !call test%run(test_getReplaced_D0_D0_D0_SK2_1, SK_"test_getReplaced_D0_D0_D0_SK2_1")
        call test%run(test_setReplaced_D0_D0_D0_SK2_1, SK_"test_setReplaced_D0_D0_D0_SK2_1")
#endif
#if     SK1_ENABLED
        !call test%run(test_getReplaced_D0_D0_D0_SK1_1, SK_"test_getReplaced_D0_D0_D0_SK1_1")
        call test%run(test_setReplaced_D0_D0_D0_SK1_1, SK_"test_setReplaced_D0_D0_D0_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_SK5_1, SK_"test_getReplaced_D1_D1_D1_SK5_1")
        call test%run(test_setReplaced_D1_D1_D1_SK5_1, SK_"test_setReplaced_D1_D1_D1_SK5_1")
#endif
#if     SK4_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_SK4_1, SK_"test_getReplaced_D1_D1_D1_SK4_1")
        call test%run(test_setReplaced_D1_D1_D1_SK4_1, SK_"test_setReplaced_D1_D1_D1_SK4_1")
#endif
#if     SK3_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_SK3_1, SK_"test_getReplaced_D1_D1_D1_SK3_1")
        call test%run(test_setReplaced_D1_D1_D1_SK3_1, SK_"test_setReplaced_D1_D1_D1_SK3_1")
#endif
#if     SK2_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_SK2_1, SK_"test_getReplaced_D1_D1_D1_SK2_1")
        call test%run(test_setReplaced_D1_D1_D1_SK2_1, SK_"test_setReplaced_D1_D1_D1_SK2_1")
#endif
#if     SK1_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_SK1_1, SK_"test_getReplaced_D1_D1_D1_SK1_1")
        call test%run(test_setReplaced_D1_D1_D1_SK1_1, SK_"test_setReplaced_D1_D1_D1_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_IK5_1, SK_"test_getReplaced_D1_D1_D1_IK5_1")
        call test%run(test_setReplaced_D1_D1_D1_IK5_1, SK_"test_setReplaced_D1_D1_D1_IK5_1")
#endif
#if     IK4_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_IK4_1, SK_"test_getReplaced_D1_D1_D1_IK4_1")
        call test%run(test_setReplaced_D1_D1_D1_IK4_1, SK_"test_setReplaced_D1_D1_D1_IK4_1")
#endif
#if     IK3_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_IK3_1, SK_"test_getReplaced_D1_D1_D1_IK3_1")
        call test%run(test_setReplaced_D1_D1_D1_IK3_1, SK_"test_setReplaced_D1_D1_D1_IK3_1")
#endif
#if     IK2_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_IK2_1, SK_"test_getReplaced_D1_D1_D1_IK2_1")
        call test%run(test_setReplaced_D1_D1_D1_IK2_1, SK_"test_setReplaced_D1_D1_D1_IK2_1")
#endif
#if     IK1_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_IK1_1, SK_"test_getReplaced_D1_D1_D1_IK1_1")
        call test%run(test_setReplaced_D1_D1_D1_IK1_1, SK_"test_setReplaced_D1_D1_D1_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_LK5_1, SK_"test_getReplaced_D1_D1_D1_LK5_1")
        call test%run(test_setReplaced_D1_D1_D1_LK5_1, SK_"test_setReplaced_D1_D1_D1_LK5_1")
#endif
#if     LK4_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_LK4_1, SK_"test_getReplaced_D1_D1_D1_LK4_1")
        call test%run(test_setReplaced_D1_D1_D1_LK4_1, SK_"test_setReplaced_D1_D1_D1_LK4_1")
#endif
#if     LK3_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_LK3_1, SK_"test_getReplaced_D1_D1_D1_LK3_1")
        call test%run(test_setReplaced_D1_D1_D1_LK3_1, SK_"test_setReplaced_D1_D1_D1_LK3_1")
#endif
#if     LK2_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_LK2_1, SK_"test_getReplaced_D1_D1_D1_LK2_1")
        call test%run(test_setReplaced_D1_D1_D1_LK2_1, SK_"test_setReplaced_D1_D1_D1_LK2_1")
#endif
#if     LK1_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_LK1_1, SK_"test_getReplaced_D1_D1_D1_LK1_1")
        call test%run(test_setReplaced_D1_D1_D1_LK1_1, SK_"test_setReplaced_D1_D1_D1_LK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_CK5_1, SK_"test_getReplaced_D1_D1_D1_CK5_1")
        call test%run(test_setReplaced_D1_D1_D1_CK5_1, SK_"test_setReplaced_D1_D1_D1_CK5_1")
#endif
#if     CK4_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_CK4_1, SK_"test_getReplaced_D1_D1_D1_CK4_1")
        call test%run(test_setReplaced_D1_D1_D1_CK4_1, SK_"test_setReplaced_D1_D1_D1_CK4_1")
#endif
#if     CK3_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_CK3_1, SK_"test_getReplaced_D1_D1_D1_CK3_1")
        call test%run(test_setReplaced_D1_D1_D1_CK3_1, SK_"test_setReplaced_D1_D1_D1_CK3_1")
#endif
#if     CK2_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_CK2_1, SK_"test_getReplaced_D1_D1_D1_CK2_1")
        call test%run(test_setReplaced_D1_D1_D1_CK2_1, SK_"test_setReplaced_D1_D1_D1_CK2_1")
#endif
#if     CK1_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_CK1_1, SK_"test_getReplaced_D1_D1_D1_CK1_1")
        call test%run(test_setReplaced_D1_D1_D1_CK1_1, SK_"test_setReplaced_D1_D1_D1_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_RK5_1, SK_"test_getReplaced_D1_D1_D1_RK5_1")
        call test%run(test_setReplaced_D1_D1_D1_RK5_1, SK_"test_setReplaced_D1_D1_D1_RK5_1")
#endif
#if     RK4_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_RK4_1, SK_"test_getReplaced_D1_D1_D1_RK4_1")
        call test%run(test_setReplaced_D1_D1_D1_RK4_1, SK_"test_setReplaced_D1_D1_D1_RK4_1")
#endif
#if     RK3_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_RK3_1, SK_"test_getReplaced_D1_D1_D1_RK3_1")
        call test%run(test_setReplaced_D1_D1_D1_RK3_1, SK_"test_setReplaced_D1_D1_D1_RK3_1")
#endif
#if     RK2_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_RK2_1, SK_"test_getReplaced_D1_D1_D1_RK2_1")
        call test%run(test_setReplaced_D1_D1_D1_RK2_1, SK_"test_setReplaced_D1_D1_D1_RK2_1")
#endif
#if     RK1_ENABLED
        !call test%run(test_getReplaced_D1_D1_D1_RK1_1, SK_"test_getReplaced_D1_D1_D1_RK1_1")
        call test%run(test_setReplaced_D1_D1_D1_RK1_1, SK_"test_setReplaced_D1_D1_D1_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_SK5_1, SK_"test_getReplaced_D1_D1_D0_SK5_1")
        call test%run(test_setReplaced_D1_D1_D0_SK5_1, SK_"test_setReplaced_D1_D1_D0_SK5_1")
#endif
#if     SK4_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_SK4_1, SK_"test_getReplaced_D1_D1_D0_SK4_1")
        call test%run(test_setReplaced_D1_D1_D0_SK4_1, SK_"test_setReplaced_D1_D1_D0_SK4_1")
#endif
#if     SK3_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_SK3_1, SK_"test_getReplaced_D1_D1_D0_SK3_1")
        call test%run(test_setReplaced_D1_D1_D0_SK3_1, SK_"test_setReplaced_D1_D1_D0_SK3_1")
#endif
#if     SK2_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_SK2_1, SK_"test_getReplaced_D1_D1_D0_SK2_1")
        call test%run(test_setReplaced_D1_D1_D0_SK2_1, SK_"test_setReplaced_D1_D1_D0_SK2_1")
#endif
#if     SK1_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_SK1_1, SK_"test_getReplaced_D1_D1_D0_SK1_1")
        call test%run(test_setReplaced_D1_D1_D0_SK1_1, SK_"test_setReplaced_D1_D1_D0_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_IK5_1, SK_"test_getReplaced_D1_D1_D0_IK5_1")
        call test%run(test_setReplaced_D1_D1_D0_IK5_1, SK_"test_setReplaced_D1_D1_D0_IK5_1")
#endif
#if     IK4_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_IK4_1, SK_"test_getReplaced_D1_D1_D0_IK4_1")
        call test%run(test_setReplaced_D1_D1_D0_IK4_1, SK_"test_setReplaced_D1_D1_D0_IK4_1")
#endif
#if     IK3_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_IK3_1, SK_"test_getReplaced_D1_D1_D0_IK3_1")
        call test%run(test_setReplaced_D1_D1_D0_IK3_1, SK_"test_setReplaced_D1_D1_D0_IK3_1")
#endif
#if     IK2_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_IK2_1, SK_"test_getReplaced_D1_D1_D0_IK2_1")
        call test%run(test_setReplaced_D1_D1_D0_IK2_1, SK_"test_setReplaced_D1_D1_D0_IK2_1")
#endif
#if     IK1_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_IK1_1, SK_"test_getReplaced_D1_D1_D0_IK1_1")
        call test%run(test_setReplaced_D1_D1_D0_IK1_1, SK_"test_setReplaced_D1_D1_D0_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_LK5_1, SK_"test_getReplaced_D1_D1_D0_LK5_1")
        call test%run(test_setReplaced_D1_D1_D0_LK5_1, SK_"test_setReplaced_D1_D1_D0_LK5_1")
#endif
#if     LK4_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_LK4_1, SK_"test_getReplaced_D1_D1_D0_LK4_1")
        call test%run(test_setReplaced_D1_D1_D0_LK4_1, SK_"test_setReplaced_D1_D1_D0_LK4_1")
#endif
#if     LK3_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_LK3_1, SK_"test_getReplaced_D1_D1_D0_LK3_1")
        call test%run(test_setReplaced_D1_D1_D0_LK3_1, SK_"test_setReplaced_D1_D1_D0_LK3_1")
#endif
#if     LK2_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_LK2_1, SK_"test_getReplaced_D1_D1_D0_LK2_1")
        call test%run(test_setReplaced_D1_D1_D0_LK2_1, SK_"test_setReplaced_D1_D1_D0_LK2_1")
#endif
#if     LK1_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_LK1_1, SK_"test_getReplaced_D1_D1_D0_LK1_1")
        call test%run(test_setReplaced_D1_D1_D0_LK1_1, SK_"test_setReplaced_D1_D1_D0_LK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_CK5_1, SK_"test_getReplaced_D1_D1_D0_CK5_1")
        call test%run(test_setReplaced_D1_D1_D0_CK5_1, SK_"test_setReplaced_D1_D1_D0_CK5_1")
#endif
#if     CK4_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_CK4_1, SK_"test_getReplaced_D1_D1_D0_CK4_1")
        call test%run(test_setReplaced_D1_D1_D0_CK4_1, SK_"test_setReplaced_D1_D1_D0_CK4_1")
#endif
#if     CK3_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_CK3_1, SK_"test_getReplaced_D1_D1_D0_CK3_1")
        call test%run(test_setReplaced_D1_D1_D0_CK3_1, SK_"test_setReplaced_D1_D1_D0_CK3_1")
#endif
#if     CK2_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_CK2_1, SK_"test_getReplaced_D1_D1_D0_CK2_1")
        call test%run(test_setReplaced_D1_D1_D0_CK2_1, SK_"test_setReplaced_D1_D1_D0_CK2_1")
#endif
#if     CK1_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_CK1_1, SK_"test_getReplaced_D1_D1_D0_CK1_1")
        call test%run(test_setReplaced_D1_D1_D0_CK1_1, SK_"test_setReplaced_D1_D1_D0_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_RK5_1, SK_"test_getReplaced_D1_D1_D0_RK5_1")
        call test%run(test_setReplaced_D1_D1_D0_RK5_1, SK_"test_setReplaced_D1_D1_D0_RK5_1")
#endif
#if     RK4_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_RK4_1, SK_"test_getReplaced_D1_D1_D0_RK4_1")
        call test%run(test_setReplaced_D1_D1_D0_RK4_1, SK_"test_setReplaced_D1_D1_D0_RK4_1")
#endif
#if     RK3_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_RK3_1, SK_"test_getReplaced_D1_D1_D0_RK3_1")
        call test%run(test_setReplaced_D1_D1_D0_RK3_1, SK_"test_setReplaced_D1_D1_D0_RK3_1")
#endif
#if     RK2_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_RK2_1, SK_"test_getReplaced_D1_D1_D0_RK2_1")
        call test%run(test_setReplaced_D1_D1_D0_RK2_1, SK_"test_setReplaced_D1_D1_D0_RK2_1")
#endif
#if     RK1_ENABLED
        !call test%run(test_getReplaced_D1_D1_D0_RK1_1, SK_"test_getReplaced_D1_D1_D0_RK1_1")
        call test%run(test_setReplaced_D1_D1_D0_RK1_1, SK_"test_setReplaced_D1_D1_D0_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_SK5_1, SK_"test_getReplaced_D1_D0_D1_SK5_1")
        call test%run(test_setReplaced_D1_D0_D1_SK5_1, SK_"test_setReplaced_D1_D0_D1_SK5_1")
#endif
#if     SK4_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_SK4_1, SK_"test_getReplaced_D1_D0_D1_SK4_1")
        call test%run(test_setReplaced_D1_D0_D1_SK4_1, SK_"test_setReplaced_D1_D0_D1_SK4_1")
#endif
#if     SK3_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_SK3_1, SK_"test_getReplaced_D1_D0_D1_SK3_1")
        call test%run(test_setReplaced_D1_D0_D1_SK3_1, SK_"test_setReplaced_D1_D0_D1_SK3_1")
#endif
#if     SK2_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_SK2_1, SK_"test_getReplaced_D1_D0_D1_SK2_1")
        call test%run(test_setReplaced_D1_D0_D1_SK2_1, SK_"test_setReplaced_D1_D0_D1_SK2_1")
#endif
#if     SK1_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_SK1_1, SK_"test_getReplaced_D1_D0_D1_SK1_1")
        call test%run(test_setReplaced_D1_D0_D1_SK1_1, SK_"test_setReplaced_D1_D0_D1_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_IK5_1, SK_"test_getReplaced_D1_D0_D1_IK5_1")
        call test%run(test_setReplaced_D1_D0_D1_IK5_1, SK_"test_setReplaced_D1_D0_D1_IK5_1")
#endif
#if     IK4_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_IK4_1, SK_"test_getReplaced_D1_D0_D1_IK4_1")
        call test%run(test_setReplaced_D1_D0_D1_IK4_1, SK_"test_setReplaced_D1_D0_D1_IK4_1")
#endif
#if     IK3_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_IK3_1, SK_"test_getReplaced_D1_D0_D1_IK3_1")
        call test%run(test_setReplaced_D1_D0_D1_IK3_1, SK_"test_setReplaced_D1_D0_D1_IK3_1")
#endif
#if     IK2_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_IK2_1, SK_"test_getReplaced_D1_D0_D1_IK2_1")
        call test%run(test_setReplaced_D1_D0_D1_IK2_1, SK_"test_setReplaced_D1_D0_D1_IK2_1")
#endif
#if     IK1_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_IK1_1, SK_"test_getReplaced_D1_D0_D1_IK1_1")
        call test%run(test_setReplaced_D1_D0_D1_IK1_1, SK_"test_setReplaced_D1_D0_D1_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_LK5_1, SK_"test_getReplaced_D1_D0_D1_LK5_1")
        call test%run(test_setReplaced_D1_D0_D1_LK5_1, SK_"test_setReplaced_D1_D0_D1_LK5_1")
#endif
#if     LK4_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_LK4_1, SK_"test_getReplaced_D1_D0_D1_LK4_1")
        call test%run(test_setReplaced_D1_D0_D1_LK4_1, SK_"test_setReplaced_D1_D0_D1_LK4_1")
#endif
#if     LK3_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_LK3_1, SK_"test_getReplaced_D1_D0_D1_LK3_1")
        call test%run(test_setReplaced_D1_D0_D1_LK3_1, SK_"test_setReplaced_D1_D0_D1_LK3_1")
#endif
#if     LK2_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_LK2_1, SK_"test_getReplaced_D1_D0_D1_LK2_1")
        call test%run(test_setReplaced_D1_D0_D1_LK2_1, SK_"test_setReplaced_D1_D0_D1_LK2_1")
#endif
#if     LK1_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_LK1_1, SK_"test_getReplaced_D1_D0_D1_LK1_1")
        call test%run(test_setReplaced_D1_D0_D1_LK1_1, SK_"test_setReplaced_D1_D0_D1_LK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_CK5_1, SK_"test_getReplaced_D1_D0_D1_CK5_1")
        call test%run(test_setReplaced_D1_D0_D1_CK5_1, SK_"test_setReplaced_D1_D0_D1_CK5_1")
#endif
#if     CK4_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_CK4_1, SK_"test_getReplaced_D1_D0_D1_CK4_1")
        call test%run(test_setReplaced_D1_D0_D1_CK4_1, SK_"test_setReplaced_D1_D0_D1_CK4_1")
#endif
#if     CK3_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_CK3_1, SK_"test_getReplaced_D1_D0_D1_CK3_1")
        call test%run(test_setReplaced_D1_D0_D1_CK3_1, SK_"test_setReplaced_D1_D0_D1_CK3_1")
#endif
#if     CK2_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_CK2_1, SK_"test_getReplaced_D1_D0_D1_CK2_1")
        call test%run(test_setReplaced_D1_D0_D1_CK2_1, SK_"test_setReplaced_D1_D0_D1_CK2_1")
#endif
#if     CK1_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_CK1_1, SK_"test_getReplaced_D1_D0_D1_CK1_1")
        call test%run(test_setReplaced_D1_D0_D1_CK1_1, SK_"test_setReplaced_D1_D0_D1_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_RK5_1, SK_"test_getReplaced_D1_D0_D1_RK5_1")
        call test%run(test_setReplaced_D1_D0_D1_RK5_1, SK_"test_setReplaced_D1_D0_D1_RK5_1")
#endif
#if     RK4_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_RK4_1, SK_"test_getReplaced_D1_D0_D1_RK4_1")
        call test%run(test_setReplaced_D1_D0_D1_RK4_1, SK_"test_setReplaced_D1_D0_D1_RK4_1")
#endif
#if     RK3_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_RK3_1, SK_"test_getReplaced_D1_D0_D1_RK3_1")
        call test%run(test_setReplaced_D1_D0_D1_RK3_1, SK_"test_setReplaced_D1_D0_D1_RK3_1")
#endif
#if     RK2_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_RK2_1, SK_"test_getReplaced_D1_D0_D1_RK2_1")
        call test%run(test_setReplaced_D1_D0_D1_RK2_1, SK_"test_setReplaced_D1_D0_D1_RK2_1")
#endif
#if     RK1_ENABLED
        !call test%run(test_getReplaced_D1_D0_D1_RK1_1, SK_"test_getReplaced_D1_D0_D1_RK1_1")
        call test%run(test_setReplaced_D1_D0_D1_RK1_1, SK_"test_setReplaced_D1_D0_D1_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_SK5_1, SK_"test_getReplaced_D1_D0_D0_SK5_1")
        call test%run(test_setReplaced_D1_D0_D0_SK5_1, SK_"test_setReplaced_D1_D0_D0_SK5_1")
#endif
#if     SK4_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_SK4_1, SK_"test_getReplaced_D1_D0_D0_SK4_1")
        call test%run(test_setReplaced_D1_D0_D0_SK4_1, SK_"test_setReplaced_D1_D0_D0_SK4_1")
#endif
#if     SK3_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_SK3_1, SK_"test_getReplaced_D1_D0_D0_SK3_1")
        call test%run(test_setReplaced_D1_D0_D0_SK3_1, SK_"test_setReplaced_D1_D0_D0_SK3_1")
#endif
#if     SK2_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_SK2_1, SK_"test_getReplaced_D1_D0_D0_SK2_1")
        call test%run(test_setReplaced_D1_D0_D0_SK2_1, SK_"test_setReplaced_D1_D0_D0_SK2_1")
#endif
#if     SK1_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_SK1_1, SK_"test_getReplaced_D1_D0_D0_SK1_1")
        call test%run(test_setReplaced_D1_D0_D0_SK1_1, SK_"test_setReplaced_D1_D0_D0_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_IK5_1, SK_"test_getReplaced_D1_D0_D0_IK5_1")
        call test%run(test_setReplaced_D1_D0_D0_IK5_1, SK_"test_setReplaced_D1_D0_D0_IK5_1")
#endif
#if     IK4_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_IK4_1, SK_"test_getReplaced_D1_D0_D0_IK4_1")
        call test%run(test_setReplaced_D1_D0_D0_IK4_1, SK_"test_setReplaced_D1_D0_D0_IK4_1")
#endif
#if     IK3_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_IK3_1, SK_"test_getReplaced_D1_D0_D0_IK3_1")
        call test%run(test_setReplaced_D1_D0_D0_IK3_1, SK_"test_setReplaced_D1_D0_D0_IK3_1")
#endif
#if     IK2_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_IK2_1, SK_"test_getReplaced_D1_D0_D0_IK2_1")
        call test%run(test_setReplaced_D1_D0_D0_IK2_1, SK_"test_setReplaced_D1_D0_D0_IK2_1")
#endif
#if     IK1_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_IK1_1, SK_"test_getReplaced_D1_D0_D0_IK1_1")
        call test%run(test_setReplaced_D1_D0_D0_IK1_1, SK_"test_setReplaced_D1_D0_D0_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_LK5_1, SK_"test_getReplaced_D1_D0_D0_LK5_1")
        call test%run(test_setReplaced_D1_D0_D0_LK5_1, SK_"test_setReplaced_D1_D0_D0_LK5_1")
#endif
#if     LK4_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_LK4_1, SK_"test_getReplaced_D1_D0_D0_LK4_1")
        call test%run(test_setReplaced_D1_D0_D0_LK4_1, SK_"test_setReplaced_D1_D0_D0_LK4_1")
#endif
#if     LK3_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_LK3_1, SK_"test_getReplaced_D1_D0_D0_LK3_1")
        call test%run(test_setReplaced_D1_D0_D0_LK3_1, SK_"test_setReplaced_D1_D0_D0_LK3_1")
#endif
#if     LK2_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_LK2_1, SK_"test_getReplaced_D1_D0_D0_LK2_1")
        call test%run(test_setReplaced_D1_D0_D0_LK2_1, SK_"test_setReplaced_D1_D0_D0_LK2_1")
#endif
#if     LK1_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_LK1_1, SK_"test_getReplaced_D1_D0_D0_LK1_1")
        call test%run(test_setReplaced_D1_D0_D0_LK1_1, SK_"test_setReplaced_D1_D0_D0_LK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_CK5_1, SK_"test_getReplaced_D1_D0_D0_CK5_1")
        call test%run(test_setReplaced_D1_D0_D0_CK5_1, SK_"test_setReplaced_D1_D0_D0_CK5_1")
#endif
#if     CK4_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_CK4_1, SK_"test_getReplaced_D1_D0_D0_CK4_1")
        call test%run(test_setReplaced_D1_D0_D0_CK4_1, SK_"test_setReplaced_D1_D0_D0_CK4_1")
#endif
#if     CK3_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_CK3_1, SK_"test_getReplaced_D1_D0_D0_CK3_1")
        call test%run(test_setReplaced_D1_D0_D0_CK3_1, SK_"test_setReplaced_D1_D0_D0_CK3_1")
#endif
#if     CK2_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_CK2_1, SK_"test_getReplaced_D1_D0_D0_CK2_1")
        call test%run(test_setReplaced_D1_D0_D0_CK2_1, SK_"test_setReplaced_D1_D0_D0_CK2_1")
#endif
#if     CK1_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_CK1_1, SK_"test_getReplaced_D1_D0_D0_CK1_1")
        call test%run(test_setReplaced_D1_D0_D0_CK1_1, SK_"test_setReplaced_D1_D0_D0_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_RK5_1, SK_"test_getReplaced_D1_D0_D0_RK5_1")
        call test%run(test_setReplaced_D1_D0_D0_RK5_1, SK_"test_setReplaced_D1_D0_D0_RK5_1")
#endif
#if     RK4_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_RK4_1, SK_"test_getReplaced_D1_D0_D0_RK4_1")
        call test%run(test_setReplaced_D1_D0_D0_RK4_1, SK_"test_setReplaced_D1_D0_D0_RK4_1")
#endif
#if     RK3_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_RK3_1, SK_"test_getReplaced_D1_D0_D0_RK3_1")
        call test%run(test_setReplaced_D1_D0_D0_RK3_1, SK_"test_setReplaced_D1_D0_D0_RK3_1")
#endif
#if     RK2_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_RK2_1, SK_"test_getReplaced_D1_D0_D0_RK2_1")
        call test%run(test_setReplaced_D1_D0_D0_RK2_1, SK_"test_setReplaced_D1_D0_D0_RK2_1")
#endif
#if     RK1_ENABLED
        !call test%run(test_getReplaced_D1_D0_D0_RK1_1, SK_"test_getReplaced_D1_D0_D0_RK1_1")
        call test%run(test_setReplaced_D1_D0_D0_RK1_1, SK_"test_setReplaced_D1_D0_D0_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%run(test_getReplaced_recurs_alloc_1, SK_"test_getReplaced_recurs_alloc_1")
        call test%run(test_getReplaced_recurs_alloc_2 , SK_"test_getReplaced_recurs_alloc_2")
        call test%run(test_getReplaced_recurs_alloc_3 , SK_"test_getReplaced_recurs_alloc_3")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getReplaced_recurs_alloc_1() result(assertion)
        implicit none
        logical(LK) :: assertion
        character(:, SK), allocatable   :: string
        string = getReplaced_recurs_alloc("  StringString ", "Str", "")
        assertion = string == "  inging "
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getReplaced_recurs_alloc_2() result(assertion)
        logical(LK)                     :: assertion
        character(:, SK), allocatable   :: string
        string = getReplaced_recurs_alloc("", "", "")
        assertion = string == ""
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getReplaced_recurs_alloc_3() result(assertion)

        implicit none
        logical(LK)                 :: assertion
        character(:, SK), allocatable   :: string
        string = getReplaced_recurs_alloc("AAA", "AA", "X")
        assertion = string == "XA"
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_arrayReplace ! LCOV_EXCL_LINE