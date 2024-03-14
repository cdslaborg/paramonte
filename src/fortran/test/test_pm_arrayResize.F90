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
!>  This module contains tests of the module [pm_arrayResize](@ref pm_arrayResize).
!>
!>  \bug
!>  \status \unresolved
!>  \source \ifort{2022}
!>  \desc
!>  There is a viscous \ifort{2022} bug where the appearance of the following `use` statements
!>  in the body of the implementation include file `test_pm_arrayRebill@routines.inc.F90` leads to various
!>  mistakes in parsing and preprocessing the contents of the include file.<br>
!>  The threshold for the maximum number of `use` statements within the entire submodule appears to be
!>  about `55`, because activating more than 55 procedures of the submodule
!>  leads to compilation failures due syntax parsing mistakes by the Intel compiler.<br>
!>  \remedy
!>  For now, all `use` statements are moved to the declaration body of the submodules.<br>
!>  This causes non-locality of the statements, but is the only solution as of today.<br>
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_arrayResize

    use pm_arrayResize
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
        module function test_setResized_D0_SK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setResized_D0_SK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setResized_D0_SK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setResized_D0_SK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setResized_D0_SK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setResized_D1_SK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setResized_D1_SK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setResized_D1_SK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setResized_D1_SK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setResized_D1_SK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setResized_D1_IK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setResized_D1_IK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setResized_D1_IK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setResized_D1_IK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setResized_D1_IK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_setResized_D1_LK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setResized_D1_LK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setResized_D1_LK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setResized_D1_LK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setResized_D1_LK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setResized_D1_CK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setResized_D1_CK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setResized_D1_CK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setResized_D1_CK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setResized_D1_CK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setResized_D1_RK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setResized_D1_RK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setResized_D1_RK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setResized_D1_RK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setResized_D1_RK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setResized_D2_SK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setResized_D2_SK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setResized_D2_SK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setResized_D2_SK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setResized_D2_SK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setResized_D2_IK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setResized_D2_IK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setResized_D2_IK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setResized_D2_IK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setResized_D2_IK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_setResized_D2_LK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setResized_D2_LK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setResized_D2_LK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setResized_D2_LK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setResized_D2_LK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setResized_D2_CK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setResized_D2_CK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setResized_D2_CK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setResized_D2_CK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setResized_D2_CK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setResized_D2_RK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setResized_D2_RK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setResized_D2_RK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setResized_D2_RK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setResized_D2_RK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setResized_D3_SK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setResized_D3_SK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setResized_D3_SK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setResized_D3_SK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setResized_D3_SK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setResized_D3_IK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setResized_D3_IK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setResized_D3_IK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setResized_D3_IK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setResized_D3_IK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_setResized_D3_LK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setResized_D3_LK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setResized_D3_LK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setResized_D3_LK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setResized_D3_LK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setResized_D3_CK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setResized_D3_CK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setResized_D3_CK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setResized_D3_CK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setResized_D3_CK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setResized_D3_RK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setResized_D3_RK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setResized_D3_RK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setResized_D3_RK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setResized_D3_RK1_1 () result(assertion); logical(LK) :: assertion; end function
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

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setResized_D0_SK5_1, SK_"test_setResized_D0_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_setResized_D0_SK4_1, SK_"test_setResized_D0_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_setResized_D0_SK3_1, SK_"test_setResized_D0_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_setResized_D0_SK2_1, SK_"test_setResized_D0_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_setResized_D0_SK1_1, SK_"test_setResized_D0_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setResized_D1_SK5_1, SK_"test_setResized_D1_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_setResized_D1_SK4_1, SK_"test_setResized_D1_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_setResized_D1_SK3_1, SK_"test_setResized_D1_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_setResized_D1_SK2_1, SK_"test_setResized_D1_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_setResized_D1_SK1_1, SK_"test_setResized_D1_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setResized_D1_IK5_1, SK_"test_setResized_D1_IK5_1")
#endif
#if     IK4_ENABLED
        call test%run(test_setResized_D1_IK4_1, SK_"test_setResized_D1_IK4_1")
#endif
#if     IK3_ENABLED
        call test%run(test_setResized_D1_IK3_1, SK_"test_setResized_D1_IK3_1")
#endif
#if     IK2_ENABLED
        call test%run(test_setResized_D1_IK2_1, SK_"test_setResized_D1_IK2_1")
#endif
#if     IK1_ENABLED
        call test%run(test_setResized_D1_IK1_1, SK_"test_setResized_D1_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setResized_D1_LK5_1, SK_"test_setResized_D1_LK5_1")
#endif
#if     LK4_ENABLED
        call test%run(test_setResized_D1_LK4_1, SK_"test_setResized_D1_LK4_1")
#endif
#if     LK3_ENABLED
        call test%run(test_setResized_D1_LK3_1, SK_"test_setResized_D1_LK3_1")
#endif
#if     LK2_ENABLED
        call test%run(test_setResized_D1_LK2_1, SK_"test_setResized_D1_LK2_1")
#endif
#if     LK1_ENABLED
        call test%run(test_setResized_D1_LK1_1, SK_"test_setResized_D1_LK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setResized_D1_CK5_1, SK_"test_setResized_D1_CK5_1")
#endif
#if     CK4_ENABLED
        call test%run(test_setResized_D1_CK4_1, SK_"test_setResized_D1_CK4_1")
#endif
#if     CK3_ENABLED
        call test%run(test_setResized_D1_CK3_1, SK_"test_setResized_D1_CK3_1")
#endif
#if     CK2_ENABLED
        call test%run(test_setResized_D1_CK2_1, SK_"test_setResized_D1_CK2_1")
#endif
#if     CK1_ENABLED
        call test%run(test_setResized_D1_CK1_1, SK_"test_setResized_D1_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setResized_D1_RK5_1, SK_"test_setResized_D1_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setResized_D1_RK4_1, SK_"test_setResized_D1_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setResized_D1_RK3_1, SK_"test_setResized_D1_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setResized_D1_RK2_1, SK_"test_setResized_D1_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setResized_D1_RK1_1, SK_"test_setResized_D1_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setResized_D2_SK5_1, SK_"test_setResized_D2_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_setResized_D2_SK4_1, SK_"test_setResized_D2_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_setResized_D2_SK3_1, SK_"test_setResized_D2_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_setResized_D2_SK2_1, SK_"test_setResized_D2_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_setResized_D2_SK1_1, SK_"test_setResized_D2_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setResized_D2_IK5_1, SK_"test_setResized_D2_IK5_1")
#endif
#if     IK4_ENABLED
        call test%run(test_setResized_D2_IK4_1, SK_"test_setResized_D2_IK4_1")
#endif
#if     IK3_ENABLED
        call test%run(test_setResized_D2_IK3_1, SK_"test_setResized_D2_IK3_1")
#endif
#if     IK2_ENABLED
        call test%run(test_setResized_D2_IK2_1, SK_"test_setResized_D2_IK2_1")
#endif
#if     IK1_ENABLED
        call test%run(test_setResized_D2_IK1_1, SK_"test_setResized_D2_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setResized_D2_LK5_1, SK_"test_setResized_D2_LK5_1")
#endif
#if     LK4_ENABLED
        call test%run(test_setResized_D2_LK4_1, SK_"test_setResized_D2_LK4_1")
#endif
#if     LK3_ENABLED
        call test%run(test_setResized_D2_LK3_1, SK_"test_setResized_D2_LK3_1")
#endif
#if     LK2_ENABLED
        call test%run(test_setResized_D2_LK2_1, SK_"test_setResized_D2_LK2_1")
#endif
#if     LK1_ENABLED
        call test%run(test_setResized_D2_LK1_1, SK_"test_setResized_D2_LK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setResized_D2_CK5_1, SK_"test_setResized_D2_CK5_1")
#endif
#if     CK4_ENABLED
        call test%run(test_setResized_D2_CK4_1, SK_"test_setResized_D2_CK4_1")
#endif
#if     CK3_ENABLED
        call test%run(test_setResized_D2_CK3_1, SK_"test_setResized_D2_CK3_1")
#endif
#if     CK2_ENABLED
        call test%run(test_setResized_D2_CK2_1, SK_"test_setResized_D2_CK2_1")
#endif
#if     CK1_ENABLED
        call test%run(test_setResized_D2_CK1_1, SK_"test_setResized_D2_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setResized_D2_RK5_1, SK_"test_setResized_D2_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setResized_D2_RK4_1, SK_"test_setResized_D2_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setResized_D2_RK3_1, SK_"test_setResized_D2_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setResized_D2_RK2_1, SK_"test_setResized_D2_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setResized_D2_RK1_1, SK_"test_setResized_D2_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setResized_D3_SK5_1, SK_"test_setResized_D3_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_setResized_D3_SK4_1, SK_"test_setResized_D3_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_setResized_D3_SK3_1, SK_"test_setResized_D3_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_setResized_D3_SK2_1, SK_"test_setResized_D3_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_setResized_D3_SK1_1, SK_"test_setResized_D3_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setResized_D3_IK5_1, SK_"test_setResized_D3_IK5_1")
#endif
#if     IK4_ENABLED
        call test%run(test_setResized_D3_IK4_1, SK_"test_setResized_D3_IK4_1")
#endif
#if     IK3_ENABLED
        call test%run(test_setResized_D3_IK3_1, SK_"test_setResized_D3_IK3_1")
#endif
#if     IK2_ENABLED
        call test%run(test_setResized_D3_IK2_1, SK_"test_setResized_D3_IK2_1")
#endif
#if     IK1_ENABLED
        call test%run(test_setResized_D3_IK1_1, SK_"test_setResized_D3_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setResized_D3_LK5_1, SK_"test_setResized_D3_LK5_1")
#endif
#if     LK4_ENABLED
        call test%run(test_setResized_D3_LK4_1, SK_"test_setResized_D3_LK4_1")
#endif
#if     LK3_ENABLED
        call test%run(test_setResized_D3_LK3_1, SK_"test_setResized_D3_LK3_1")
#endif
#if     LK2_ENABLED
        call test%run(test_setResized_D3_LK2_1, SK_"test_setResized_D3_LK2_1")
#endif
#if     LK1_ENABLED
        call test%run(test_setResized_D3_LK1_1, SK_"test_setResized_D3_LK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setResized_D3_CK5_1, SK_"test_setResized_D3_CK5_1")
#endif
#if     CK4_ENABLED
        call test%run(test_setResized_D3_CK4_1, SK_"test_setResized_D3_CK4_1")
#endif
#if     CK3_ENABLED
        call test%run(test_setResized_D3_CK3_1, SK_"test_setResized_D3_CK3_1")
#endif
#if     CK2_ENABLED
        call test%run(test_setResized_D3_CK2_1, SK_"test_setResized_D3_CK2_1")
#endif
#if     CK1_ENABLED
        call test%run(test_setResized_D3_CK1_1, SK_"test_setResized_D3_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setResized_D3_RK5_1, SK_"test_setResized_D3_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setResized_D3_RK4_1, SK_"test_setResized_D3_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setResized_D3_RK3_1, SK_"test_setResized_D3_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setResized_D3_RK2_1, SK_"test_setResized_D3_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setResized_D3_RK1_1, SK_"test_setResized_D3_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_arrayResize