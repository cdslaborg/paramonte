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
!>  This module contains tests of the module [pm_arraySort](@ref pm_arraySort).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_arraySort

    use pm_arraySort
    use pm_test, only: test_type, LK
    use pm_kind, only: IK, RK
    use pm_kind, only: LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_isAscending_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_isAscending_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_isAscending_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_isAscending_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_isAscending_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_isAscending_D1_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_isAscending_D1_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_isAscending_D1_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_isAscending_D1_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_isAscending_D1_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_isAscending_D1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_isAscending_D1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_isAscending_D1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_isAscending_D1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_isAscending_D1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_isAscending_D1_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_isAscending_D1_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_isAscending_D1_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_isAscending_D1_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_isAscending_D1_LK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_isAscending_D1_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_isAscending_D1_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_isAscending_D1_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_isAscending_D1_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_isAscending_D1_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_isAscending_D1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_isAscending_D1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_isAscending_D1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_isAscending_D1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_isAscending_D1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_isAscending_D1_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_isAscending_D1_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_isAscending_D1_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_isAscending_D1_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_isAscending_D1_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_isDescending_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_isDescending_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_isDescending_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_isDescending_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_isDescending_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_isDescending_D1_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_isDescending_D1_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_isDescending_D1_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_isDescending_D1_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_isDescending_D1_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_isDescending_D1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_isDescending_D1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_isDescending_D1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_isDescending_D1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_isDescending_D1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_isDescending_D1_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_isDescending_D1_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_isDescending_D1_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_isDescending_D1_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_isDescending_D1_LK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_isDescending_D1_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_isDescending_D1_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_isDescending_D1_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_isDescending_D1_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_isDescending_D1_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_isDescending_D1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_isDescending_D1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_isDescending_D1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_isDescending_D1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_isDescending_D1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_isDescending_D1_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_isDescending_D1_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_isDescending_D1_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_isDescending_D1_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_isDescending_D1_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_isSorted_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_isSorted_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_isSorted_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_isSorted_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_isSorted_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_isSorted_D1_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_isSorted_D1_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_isSorted_D1_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_isSorted_D1_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_isSorted_D1_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_isSorted_D1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_isSorted_D1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_isSorted_D1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_isSorted_D1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_isSorted_D1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_isSorted_D1_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_isSorted_D1_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_isSorted_D1_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_isSorted_D1_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_isSorted_D1_LK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_isSorted_D1_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_isSorted_D1_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_isSorted_D1_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_isSorted_D1_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_isSorted_D1_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_isSorted_D1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_isSorted_D1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_isSorted_D1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_isSorted_D1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_isSorted_D1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_isSorted_D1_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_isSorted_D1_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_isSorted_D1_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_isSorted_D1_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_isSorted_D1_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedIndDef_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedIndDef_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedIndDef_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedIndDef_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedIndDef_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedIndDef_D1_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedIndDef_D1_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedIndDef_D1_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedIndDef_D1_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedIndDef_D1_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setSortedIndDef_D1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setSortedIndDef_D1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setSortedIndDef_D1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setSortedIndDef_D1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setSortedIndDef_D1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_setSortedIndDef_D1_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setSortedIndDef_D1_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setSortedIndDef_D1_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setSortedIndDef_D1_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setSortedIndDef_D1_LK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setSortedIndDef_D1_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setSortedIndDef_D1_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setSortedIndDef_D1_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setSortedIndDef_D1_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setSortedIndDef_D1_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setSortedIndDef_D1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setSortedIndDef_D1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setSortedIndDef_D1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setSortedIndDef_D1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setSortedIndDef_D1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedIndDef_D1_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedIndDef_D1_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedIndDef_D1_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedIndDef_D1_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedIndDef_D1_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrQsorti_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrQsorti_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrQsorti_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrQsorti_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrQsorti_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrQsorti_D1_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrQsorti_D1_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrQsorti_D1_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrQsorti_D1_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrQsorti_D1_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setSortedArrQsorti_D1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setSortedArrQsorti_D1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setSortedArrQsorti_D1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setSortedArrQsorti_D1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setSortedArrQsorti_D1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_setSortedArrQsorti_D1_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setSortedArrQsorti_D1_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setSortedArrQsorti_D1_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setSortedArrQsorti_D1_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setSortedArrQsorti_D1_LK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setSortedArrQsorti_D1_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setSortedArrQsorti_D1_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setSortedArrQsorti_D1_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setSortedArrQsorti_D1_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setSortedArrQsorti_D1_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setSortedArrQsorti_D1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setSortedArrQsorti_D1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setSortedArrQsorti_D1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setSortedArrQsorti_D1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setSortedArrQsorti_D1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrQsorti_D1_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrQsorti_D1_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrQsorti_D1_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrQsorti_D1_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrQsorti_D1_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrQsortr_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrQsortr_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrQsortr_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrQsortr_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrQsortr_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrQsortr_D1_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrQsortr_D1_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrQsortr_D1_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrQsortr_D1_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrQsortr_D1_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setSortedArrQsortr_D1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setSortedArrQsortr_D1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setSortedArrQsortr_D1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setSortedArrQsortr_D1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setSortedArrQsortr_D1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_setSortedArrQsortr_D1_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setSortedArrQsortr_D1_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setSortedArrQsortr_D1_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setSortedArrQsortr_D1_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setSortedArrQsortr_D1_LK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setSortedArrQsortr_D1_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setSortedArrQsortr_D1_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setSortedArrQsortr_D1_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setSortedArrQsortr_D1_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setSortedArrQsortr_D1_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setSortedArrQsortr_D1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setSortedArrQsortr_D1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setSortedArrQsortr_D1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setSortedArrQsortr_D1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setSortedArrQsortr_D1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrQsortr_D1_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrQsortr_D1_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrQsortr_D1_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrQsortr_D1_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrQsortr_D1_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrQsortrdp_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrQsortrdp_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrQsortrdp_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrQsortrdp_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrQsortrdp_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrQsortrdp_D1_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrQsortrdp_D1_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrQsortrdp_D1_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrQsortrdp_D1_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrQsortrdp_D1_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setSortedArrQsortrdp_D1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setSortedArrQsortrdp_D1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setSortedArrQsortrdp_D1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setSortedArrQsortrdp_D1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setSortedArrQsortrdp_D1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_setSortedArrQsortrdp_D1_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setSortedArrQsortrdp_D1_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setSortedArrQsortrdp_D1_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setSortedArrQsortrdp_D1_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setSortedArrQsortrdp_D1_LK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setSortedArrQsortrdp_D1_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setSortedArrQsortrdp_D1_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setSortedArrQsortrdp_D1_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setSortedArrQsortrdp_D1_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setSortedArrQsortrdp_D1_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setSortedArrQsortrdp_D1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setSortedArrQsortrdp_D1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setSortedArrQsortrdp_D1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setSortedArrQsortrdp_D1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setSortedArrQsortrdp_D1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrQsortrdp_D1_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrQsortrdp_D1_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrQsortrdp_D1_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrQsortrdp_D1_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrQsortrdp_D1_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrBubble_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrBubble_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrBubble_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrBubble_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrBubble_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrBubble_D1_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrBubble_D1_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrBubble_D1_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrBubble_D1_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrBubble_D1_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setSortedArrBubble_D1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setSortedArrBubble_D1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setSortedArrBubble_D1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setSortedArrBubble_D1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setSortedArrBubble_D1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_setSortedArrBubble_D1_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setSortedArrBubble_D1_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setSortedArrBubble_D1_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setSortedArrBubble_D1_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setSortedArrBubble_D1_LK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setSortedArrBubble_D1_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setSortedArrBubble_D1_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setSortedArrBubble_D1_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setSortedArrBubble_D1_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setSortedArrBubble_D1_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setSortedArrBubble_D1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setSortedArrBubble_D1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setSortedArrBubble_D1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setSortedArrBubble_D1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setSortedArrBubble_D1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrBubble_D1_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrBubble_D1_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrBubble_D1_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrBubble_D1_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrBubble_D1_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrHeapi_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrHeapi_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrHeapi_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrHeapi_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrHeapi_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrHeapi_D1_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrHeapi_D1_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrHeapi_D1_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrHeapi_D1_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrHeapi_D1_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setSortedArrHeapi_D1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setSortedArrHeapi_D1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setSortedArrHeapi_D1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setSortedArrHeapi_D1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setSortedArrHeapi_D1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_setSortedArrHeapi_D1_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setSortedArrHeapi_D1_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setSortedArrHeapi_D1_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setSortedArrHeapi_D1_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setSortedArrHeapi_D1_LK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setSortedArrHeapi_D1_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setSortedArrHeapi_D1_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setSortedArrHeapi_D1_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setSortedArrHeapi_D1_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setSortedArrHeapi_D1_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setSortedArrHeapi_D1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setSortedArrHeapi_D1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setSortedArrHeapi_D1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setSortedArrHeapi_D1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setSortedArrHeapi_D1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrHeapi_D1_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrHeapi_D1_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrHeapi_D1_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrHeapi_D1_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrHeapi_D1_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrHeapr_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrHeapr_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrHeapr_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrHeapr_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrHeapr_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrHeapr_D1_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrHeapr_D1_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrHeapr_D1_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrHeapr_D1_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrHeapr_D1_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setSortedArrHeapr_D1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setSortedArrHeapr_D1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setSortedArrHeapr_D1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setSortedArrHeapr_D1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setSortedArrHeapr_D1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_setSortedArrHeapr_D1_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setSortedArrHeapr_D1_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setSortedArrHeapr_D1_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setSortedArrHeapr_D1_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setSortedArrHeapr_D1_LK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setSortedArrHeapr_D1_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setSortedArrHeapr_D1_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setSortedArrHeapr_D1_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setSortedArrHeapr_D1_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setSortedArrHeapr_D1_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setSortedArrHeapr_D1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setSortedArrHeapr_D1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setSortedArrHeapr_D1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setSortedArrHeapr_D1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setSortedArrHeapr_D1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrHeapr_D1_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrHeapr_D1_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrHeapr_D1_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrHeapr_D1_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrHeapr_D1_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrInsertionl_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrInsertionl_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrInsertionl_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrInsertionl_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrInsertionl_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrInsertionl_D1_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrInsertionl_D1_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrInsertionl_D1_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrInsertionl_D1_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrInsertionl_D1_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setSortedArrInsertionl_D1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setSortedArrInsertionl_D1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setSortedArrInsertionl_D1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setSortedArrInsertionl_D1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setSortedArrInsertionl_D1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_setSortedArrInsertionl_D1_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setSortedArrInsertionl_D1_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setSortedArrInsertionl_D1_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setSortedArrInsertionl_D1_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setSortedArrInsertionl_D1_LK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setSortedArrInsertionl_D1_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setSortedArrInsertionl_D1_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setSortedArrInsertionl_D1_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setSortedArrInsertionl_D1_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setSortedArrInsertionl_D1_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setSortedArrInsertionl_D1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setSortedArrInsertionl_D1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setSortedArrInsertionl_D1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setSortedArrInsertionl_D1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setSortedArrInsertionl_D1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrInsertionl_D1_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrInsertionl_D1_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrInsertionl_D1_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrInsertionl_D1_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrInsertionl_D1_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrInsertionb_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrInsertionb_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrInsertionb_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrInsertionb_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrInsertionb_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrInsertionb_D1_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrInsertionb_D1_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrInsertionb_D1_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrInsertionb_D1_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrInsertionb_D1_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setSortedArrInsertionb_D1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setSortedArrInsertionb_D1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setSortedArrInsertionb_D1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setSortedArrInsertionb_D1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setSortedArrInsertionb_D1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_setSortedArrInsertionb_D1_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setSortedArrInsertionb_D1_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setSortedArrInsertionb_D1_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setSortedArrInsertionb_D1_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setSortedArrInsertionb_D1_LK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setSortedArrInsertionb_D1_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setSortedArrInsertionb_D1_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setSortedArrInsertionb_D1_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setSortedArrInsertionb_D1_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setSortedArrInsertionb_D1_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setSortedArrInsertionb_D1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setSortedArrInsertionb_D1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setSortedArrInsertionb_D1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setSortedArrInsertionb_D1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setSortedArrInsertionb_D1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrInsertionb_D1_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrInsertionb_D1_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrInsertionb_D1_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrInsertionb_D1_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrInsertionb_D1_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrMerger_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrMerger_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrMerger_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrMerger_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrMerger_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrMerger_D1_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrMerger_D1_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrMerger_D1_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrMerger_D1_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrMerger_D1_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setSortedArrMerger_D1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setSortedArrMerger_D1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setSortedArrMerger_D1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setSortedArrMerger_D1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setSortedArrMerger_D1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_setSortedArrMerger_D1_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setSortedArrMerger_D1_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setSortedArrMerger_D1_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setSortedArrMerger_D1_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setSortedArrMerger_D1_LK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setSortedArrMerger_D1_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setSortedArrMerger_D1_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setSortedArrMerger_D1_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setSortedArrMerger_D1_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setSortedArrMerger_D1_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setSortedArrMerger_D1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setSortedArrMerger_D1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setSortedArrMerger_D1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setSortedArrMerger_D1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setSortedArrMerger_D1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrMerger_D1_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrMerger_D1_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrMerger_D1_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrMerger_D1_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrMerger_D1_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrSelection_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrSelection_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrSelection_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrSelection_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrSelection_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrSelection_D1_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrSelection_D1_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrSelection_D1_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrSelection_D1_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrSelection_D1_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setSortedArrSelection_D1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setSortedArrSelection_D1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setSortedArrSelection_D1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setSortedArrSelection_D1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setSortedArrSelection_D1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_setSortedArrSelection_D1_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setSortedArrSelection_D1_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setSortedArrSelection_D1_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setSortedArrSelection_D1_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setSortedArrSelection_D1_LK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setSortedArrSelection_D1_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setSortedArrSelection_D1_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setSortedArrSelection_D1_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setSortedArrSelection_D1_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setSortedArrSelection_D1_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setSortedArrSelection_D1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setSortedArrSelection_D1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setSortedArrSelection_D1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setSortedArrSelection_D1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setSortedArrSelection_D1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrSelection_D1_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrSelection_D1_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrSelection_D1_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrSelection_D1_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrSelection_D1_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrShell_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrShell_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrShell_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrShell_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrShell_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrShell_D1_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrShell_D1_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrShell_D1_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrShell_D1_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrShell_D1_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setSortedArrShell_D1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setSortedArrShell_D1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setSortedArrShell_D1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setSortedArrShell_D1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setSortedArrShell_D1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_setSortedArrShell_D1_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setSortedArrShell_D1_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setSortedArrShell_D1_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setSortedArrShell_D1_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setSortedArrShell_D1_LK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setSortedArrShell_D1_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setSortedArrShell_D1_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setSortedArrShell_D1_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setSortedArrShell_D1_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setSortedArrShell_D1_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setSortedArrShell_D1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setSortedArrShell_D1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setSortedArrShell_D1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setSortedArrShell_D1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setSortedArrShell_D1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setSortedArrShell_D1_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setSortedArrShell_D1_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setSortedArrShell_D1_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setSortedArrShell_D1_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setSortedArrShell_D1_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()
        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_isAscending_D0_SK5, SK_"test_isAscending_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_isAscending_D0_SK4, SK_"test_isAscending_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_isAscending_D0_SK3, SK_"test_isAscending_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_isAscending_D0_SK2, SK_"test_isAscending_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_isAscending_D0_SK1, SK_"test_isAscending_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_isAscending_D1_SK5, SK_"test_isAscending_D1_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_isAscending_D1_SK4, SK_"test_isAscending_D1_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_isAscending_D1_SK3, SK_"test_isAscending_D1_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_isAscending_D1_SK2, SK_"test_isAscending_D1_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_isAscending_D1_SK1, SK_"test_isAscending_D1_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_isAscending_D1_IK5, SK_"test_isAscending_D1_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_isAscending_D1_IK4, SK_"test_isAscending_D1_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_isAscending_D1_IK3, SK_"test_isAscending_D1_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_isAscending_D1_IK2, SK_"test_isAscending_D1_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_isAscending_D1_IK1, SK_"test_isAscending_D1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_isAscending_D1_LK5, SK_"test_isAscending_D1_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_isAscending_D1_LK4, SK_"test_isAscending_D1_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_isAscending_D1_LK3, SK_"test_isAscending_D1_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_isAscending_D1_LK2, SK_"test_isAscending_D1_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_isAscending_D1_LK1, SK_"test_isAscending_D1_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_isAscending_D1_CK5, SK_"test_isAscending_D1_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_isAscending_D1_CK4, SK_"test_isAscending_D1_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_isAscending_D1_CK3, SK_"test_isAscending_D1_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_isAscending_D1_CK2, SK_"test_isAscending_D1_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_isAscending_D1_CK1, SK_"test_isAscending_D1_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_isAscending_D1_RK5, SK_"test_isAscending_D1_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_isAscending_D1_RK4, SK_"test_isAscending_D1_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_isAscending_D1_RK3, SK_"test_isAscending_D1_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_isAscending_D1_RK2, SK_"test_isAscending_D1_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_isAscending_D1_RK1, SK_"test_isAscending_D1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
#if     SK5_ENABLED
        call test%run(test_isAscending_D1_PSSK5, SK_"test_isAscending_D1_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_isAscending_D1_PSSK4, SK_"test_isAscending_D1_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_isAscending_D1_PSSK3, SK_"test_isAscending_D1_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_isAscending_D1_PSSK2, SK_"test_isAscending_D1_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_isAscending_D1_PSSK1, SK_"test_isAscending_D1_PSSK1")
#endif
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_isDescending_D0_SK5, SK_"test_isDescending_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_isDescending_D0_SK4, SK_"test_isDescending_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_isDescending_D0_SK3, SK_"test_isDescending_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_isDescending_D0_SK2, SK_"test_isDescending_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_isDescending_D0_SK1, SK_"test_isDescending_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_isDescending_D1_SK5, SK_"test_isDescending_D1_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_isDescending_D1_SK4, SK_"test_isDescending_D1_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_isDescending_D1_SK3, SK_"test_isDescending_D1_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_isDescending_D1_SK2, SK_"test_isDescending_D1_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_isDescending_D1_SK1, SK_"test_isDescending_D1_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_isDescending_D1_IK5, SK_"test_isDescending_D1_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_isDescending_D1_IK4, SK_"test_isDescending_D1_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_isDescending_D1_IK3, SK_"test_isDescending_D1_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_isDescending_D1_IK2, SK_"test_isDescending_D1_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_isDescending_D1_IK1, SK_"test_isDescending_D1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_isDescending_D1_LK5, SK_"test_isDescending_D1_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_isDescending_D1_LK4, SK_"test_isDescending_D1_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_isDescending_D1_LK3, SK_"test_isDescending_D1_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_isDescending_D1_LK2, SK_"test_isDescending_D1_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_isDescending_D1_LK1, SK_"test_isDescending_D1_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_isDescending_D1_CK5, SK_"test_isDescending_D1_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_isDescending_D1_CK4, SK_"test_isDescending_D1_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_isDescending_D1_CK3, SK_"test_isDescending_D1_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_isDescending_D1_CK2, SK_"test_isDescending_D1_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_isDescending_D1_CK1, SK_"test_isDescending_D1_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_isDescending_D1_RK5, SK_"test_isDescending_D1_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_isDescending_D1_RK4, SK_"test_isDescending_D1_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_isDescending_D1_RK3, SK_"test_isDescending_D1_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_isDescending_D1_RK2, SK_"test_isDescending_D1_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_isDescending_D1_RK1, SK_"test_isDescending_D1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
#if     SK5_ENABLED
        call test%run(test_isDescending_D1_PSSK5, SK_"test_isDescending_D1_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_isDescending_D1_PSSK4, SK_"test_isDescending_D1_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_isDescending_D1_PSSK3, SK_"test_isDescending_D1_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_isDescending_D1_PSSK2, SK_"test_isDescending_D1_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_isDescending_D1_PSSK1, SK_"test_isDescending_D1_PSSK1")
#endif
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_isSorted_D0_SK5, SK_"test_isSorted_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_isSorted_D0_SK4, SK_"test_isSorted_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_isSorted_D0_SK3, SK_"test_isSorted_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_isSorted_D0_SK2, SK_"test_isSorted_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_isSorted_D0_SK1, SK_"test_isSorted_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_isSorted_D1_SK5, SK_"test_isSorted_D1_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_isSorted_D1_SK4, SK_"test_isSorted_D1_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_isSorted_D1_SK3, SK_"test_isSorted_D1_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_isSorted_D1_SK2, SK_"test_isSorted_D1_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_isSorted_D1_SK1, SK_"test_isSorted_D1_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_isSorted_D1_IK5, SK_"test_isSorted_D1_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_isSorted_D1_IK4, SK_"test_isSorted_D1_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_isSorted_D1_IK3, SK_"test_isSorted_D1_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_isSorted_D1_IK2, SK_"test_isSorted_D1_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_isSorted_D1_IK1, SK_"test_isSorted_D1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_isSorted_D1_LK5, SK_"test_isSorted_D1_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_isSorted_D1_LK4, SK_"test_isSorted_D1_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_isSorted_D1_LK3, SK_"test_isSorted_D1_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_isSorted_D1_LK2, SK_"test_isSorted_D1_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_isSorted_D1_LK1, SK_"test_isSorted_D1_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_isSorted_D1_CK5, SK_"test_isSorted_D1_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_isSorted_D1_CK4, SK_"test_isSorted_D1_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_isSorted_D1_CK3, SK_"test_isSorted_D1_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_isSorted_D1_CK2, SK_"test_isSorted_D1_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_isSorted_D1_CK1, SK_"test_isSorted_D1_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_isSorted_D1_RK5, SK_"test_isSorted_D1_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_isSorted_D1_RK4, SK_"test_isSorted_D1_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_isSorted_D1_RK3, SK_"test_isSorted_D1_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_isSorted_D1_RK2, SK_"test_isSorted_D1_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_isSorted_D1_RK1, SK_"test_isSorted_D1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
#if     SK5_ENABLED
        call test%run(test_isSorted_D1_PSSK5, SK_"test_isSorted_D1_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_isSorted_D1_PSSK4, SK_"test_isSorted_D1_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_isSorted_D1_PSSK3, SK_"test_isSorted_D1_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_isSorted_D1_PSSK2, SK_"test_isSorted_D1_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_isSorted_D1_PSSK1, SK_"test_isSorted_D1_PSSK1")
#endif
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedIndDef_D0_SK5, SK_"test_setSortedIndDef_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedIndDef_D0_SK4, SK_"test_setSortedIndDef_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedIndDef_D0_SK3, SK_"test_setSortedIndDef_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedIndDef_D0_SK2, SK_"test_setSortedIndDef_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedIndDef_D0_SK1, SK_"test_setSortedIndDef_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedIndDef_D1_SK5, SK_"test_setSortedIndDef_D1_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedIndDef_D1_SK4, SK_"test_setSortedIndDef_D1_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedIndDef_D1_SK3, SK_"test_setSortedIndDef_D1_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedIndDef_D1_SK2, SK_"test_setSortedIndDef_D1_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedIndDef_D1_SK1, SK_"test_setSortedIndDef_D1_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setSortedIndDef_D1_IK5, SK_"test_setSortedIndDef_D1_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_setSortedIndDef_D1_IK4, SK_"test_setSortedIndDef_D1_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_setSortedIndDef_D1_IK3, SK_"test_setSortedIndDef_D1_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_setSortedIndDef_D1_IK2, SK_"test_setSortedIndDef_D1_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_setSortedIndDef_D1_IK1, SK_"test_setSortedIndDef_D1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setSortedIndDef_D1_LK5, SK_"test_setSortedIndDef_D1_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_setSortedIndDef_D1_LK4, SK_"test_setSortedIndDef_D1_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_setSortedIndDef_D1_LK3, SK_"test_setSortedIndDef_D1_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_setSortedIndDef_D1_LK2, SK_"test_setSortedIndDef_D1_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_setSortedIndDef_D1_LK1, SK_"test_setSortedIndDef_D1_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setSortedIndDef_D1_CK5, SK_"test_setSortedIndDef_D1_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setSortedIndDef_D1_CK4, SK_"test_setSortedIndDef_D1_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setSortedIndDef_D1_CK3, SK_"test_setSortedIndDef_D1_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setSortedIndDef_D1_CK2, SK_"test_setSortedIndDef_D1_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setSortedIndDef_D1_CK1, SK_"test_setSortedIndDef_D1_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setSortedIndDef_D1_RK5, SK_"test_setSortedIndDef_D1_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setSortedIndDef_D1_RK4, SK_"test_setSortedIndDef_D1_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setSortedIndDef_D1_RK3, SK_"test_setSortedIndDef_D1_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setSortedIndDef_D1_RK2, SK_"test_setSortedIndDef_D1_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setSortedIndDef_D1_RK1, SK_"test_setSortedIndDef_D1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
#if     SK5_ENABLED
        call test%run(test_setSortedIndDef_D1_PSSK5, SK_"test_setSortedIndDef_D1_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedIndDef_D1_PSSK4, SK_"test_setSortedIndDef_D1_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedIndDef_D1_PSSK3, SK_"test_setSortedIndDef_D1_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedIndDef_D1_PSSK2, SK_"test_setSortedIndDef_D1_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedIndDef_D1_PSSK1, SK_"test_setSortedIndDef_D1_PSSK1")
#endif
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrQsorti_D0_SK5, SK_"test_setSortedArrQsorti_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrQsorti_D0_SK4, SK_"test_setSortedArrQsorti_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrQsorti_D0_SK3, SK_"test_setSortedArrQsorti_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrQsorti_D0_SK2, SK_"test_setSortedArrQsorti_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrQsorti_D0_SK1, SK_"test_setSortedArrQsorti_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrQsorti_D1_SK5, SK_"test_setSortedArrQsorti_D1_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrQsorti_D1_SK4, SK_"test_setSortedArrQsorti_D1_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrQsorti_D1_SK3, SK_"test_setSortedArrQsorti_D1_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrQsorti_D1_SK2, SK_"test_setSortedArrQsorti_D1_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrQsorti_D1_SK1, SK_"test_setSortedArrQsorti_D1_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setSortedArrQsorti_D1_IK5, SK_"test_setSortedArrQsorti_D1_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_setSortedArrQsorti_D1_IK4, SK_"test_setSortedArrQsorti_D1_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_setSortedArrQsorti_D1_IK3, SK_"test_setSortedArrQsorti_D1_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_setSortedArrQsorti_D1_IK2, SK_"test_setSortedArrQsorti_D1_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_setSortedArrQsorti_D1_IK1, SK_"test_setSortedArrQsorti_D1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setSortedArrQsorti_D1_LK5, SK_"test_setSortedArrQsorti_D1_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_setSortedArrQsorti_D1_LK4, SK_"test_setSortedArrQsorti_D1_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_setSortedArrQsorti_D1_LK3, SK_"test_setSortedArrQsorti_D1_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_setSortedArrQsorti_D1_LK2, SK_"test_setSortedArrQsorti_D1_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_setSortedArrQsorti_D1_LK1, SK_"test_setSortedArrQsorti_D1_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setSortedArrQsorti_D1_CK5, SK_"test_setSortedArrQsorti_D1_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setSortedArrQsorti_D1_CK4, SK_"test_setSortedArrQsorti_D1_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setSortedArrQsorti_D1_CK3, SK_"test_setSortedArrQsorti_D1_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setSortedArrQsorti_D1_CK2, SK_"test_setSortedArrQsorti_D1_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setSortedArrQsorti_D1_CK1, SK_"test_setSortedArrQsorti_D1_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setSortedArrQsorti_D1_RK5, SK_"test_setSortedArrQsorti_D1_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setSortedArrQsorti_D1_RK4, SK_"test_setSortedArrQsorti_D1_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setSortedArrQsorti_D1_RK3, SK_"test_setSortedArrQsorti_D1_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setSortedArrQsorti_D1_RK2, SK_"test_setSortedArrQsorti_D1_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setSortedArrQsorti_D1_RK1, SK_"test_setSortedArrQsorti_D1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
#if     SK5_ENABLED
        call test%run(test_setSortedArrQsorti_D1_PSSK5, SK_"test_setSortedArrQsorti_D1_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrQsorti_D1_PSSK4, SK_"test_setSortedArrQsorti_D1_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrQsorti_D1_PSSK3, SK_"test_setSortedArrQsorti_D1_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrQsorti_D1_PSSK2, SK_"test_setSortedArrQsorti_D1_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrQsorti_D1_PSSK1, SK_"test_setSortedArrQsorti_D1_PSSK1")
#endif
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrQsortr_D0_SK5, SK_"test_setSortedArrQsortr_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrQsortr_D0_SK4, SK_"test_setSortedArrQsortr_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrQsortr_D0_SK3, SK_"test_setSortedArrQsortr_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrQsortr_D0_SK2, SK_"test_setSortedArrQsortr_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrQsortr_D0_SK1, SK_"test_setSortedArrQsortr_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrQsortr_D1_SK5, SK_"test_setSortedArrQsortr_D1_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrQsortr_D1_SK4, SK_"test_setSortedArrQsortr_D1_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrQsortr_D1_SK3, SK_"test_setSortedArrQsortr_D1_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrQsortr_D1_SK2, SK_"test_setSortedArrQsortr_D1_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrQsortr_D1_SK1, SK_"test_setSortedArrQsortr_D1_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setSortedArrQsortr_D1_IK5, SK_"test_setSortedArrQsortr_D1_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_setSortedArrQsortr_D1_IK4, SK_"test_setSortedArrQsortr_D1_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_setSortedArrQsortr_D1_IK3, SK_"test_setSortedArrQsortr_D1_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_setSortedArrQsortr_D1_IK2, SK_"test_setSortedArrQsortr_D1_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_setSortedArrQsortr_D1_IK1, SK_"test_setSortedArrQsortr_D1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setSortedArrQsortr_D1_LK5, SK_"test_setSortedArrQsortr_D1_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_setSortedArrQsortr_D1_LK4, SK_"test_setSortedArrQsortr_D1_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_setSortedArrQsortr_D1_LK3, SK_"test_setSortedArrQsortr_D1_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_setSortedArrQsortr_D1_LK2, SK_"test_setSortedArrQsortr_D1_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_setSortedArrQsortr_D1_LK1, SK_"test_setSortedArrQsortr_D1_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setSortedArrQsortr_D1_CK5, SK_"test_setSortedArrQsortr_D1_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setSortedArrQsortr_D1_CK4, SK_"test_setSortedArrQsortr_D1_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setSortedArrQsortr_D1_CK3, SK_"test_setSortedArrQsortr_D1_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setSortedArrQsortr_D1_CK2, SK_"test_setSortedArrQsortr_D1_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setSortedArrQsortr_D1_CK1, SK_"test_setSortedArrQsortr_D1_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setSortedArrQsortr_D1_RK5, SK_"test_setSortedArrQsortr_D1_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setSortedArrQsortr_D1_RK4, SK_"test_setSortedArrQsortr_D1_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setSortedArrQsortr_D1_RK3, SK_"test_setSortedArrQsortr_D1_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setSortedArrQsortr_D1_RK2, SK_"test_setSortedArrQsortr_D1_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setSortedArrQsortr_D1_RK1, SK_"test_setSortedArrQsortr_D1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
#if     SK5_ENABLED
        call test%run(test_setSortedArrQsortr_D1_PSSK5, SK_"test_setSortedArrQsortr_D1_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrQsortr_D1_PSSK4, SK_"test_setSortedArrQsortr_D1_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrQsortr_D1_PSSK3, SK_"test_setSortedArrQsortr_D1_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrQsortr_D1_PSSK2, SK_"test_setSortedArrQsortr_D1_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrQsortr_D1_PSSK1, SK_"test_setSortedArrQsortr_D1_PSSK1")
#endif
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrQsortrdp_D0_SK5, SK_"test_setSortedArrQsortrdp_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrQsortrdp_D0_SK4, SK_"test_setSortedArrQsortrdp_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrQsortrdp_D0_SK3, SK_"test_setSortedArrQsortrdp_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrQsortrdp_D0_SK2, SK_"test_setSortedArrQsortrdp_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrQsortrdp_D0_SK1, SK_"test_setSortedArrQsortrdp_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_SK5, SK_"test_setSortedArrQsortrdp_D1_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_SK4, SK_"test_setSortedArrQsortrdp_D1_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_SK3, SK_"test_setSortedArrQsortrdp_D1_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_SK2, SK_"test_setSortedArrQsortrdp_D1_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_SK1, SK_"test_setSortedArrQsortrdp_D1_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_IK5, SK_"test_setSortedArrQsortrdp_D1_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_IK4, SK_"test_setSortedArrQsortrdp_D1_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_IK3, SK_"test_setSortedArrQsortrdp_D1_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_IK2, SK_"test_setSortedArrQsortrdp_D1_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_IK1, SK_"test_setSortedArrQsortrdp_D1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_LK5, SK_"test_setSortedArrQsortrdp_D1_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_LK4, SK_"test_setSortedArrQsortrdp_D1_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_LK3, SK_"test_setSortedArrQsortrdp_D1_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_LK2, SK_"test_setSortedArrQsortrdp_D1_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_LK1, SK_"test_setSortedArrQsortrdp_D1_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_CK5, SK_"test_setSortedArrQsortrdp_D1_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_CK4, SK_"test_setSortedArrQsortrdp_D1_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_CK3, SK_"test_setSortedArrQsortrdp_D1_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_CK2, SK_"test_setSortedArrQsortrdp_D1_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_CK1, SK_"test_setSortedArrQsortrdp_D1_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_RK5, SK_"test_setSortedArrQsortrdp_D1_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_RK4, SK_"test_setSortedArrQsortrdp_D1_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_RK3, SK_"test_setSortedArrQsortrdp_D1_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_RK2, SK_"test_setSortedArrQsortrdp_D1_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_RK1, SK_"test_setSortedArrQsortrdp_D1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
#if     SK5_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_PSSK5, SK_"test_setSortedArrQsortrdp_D1_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_PSSK4, SK_"test_setSortedArrQsortrdp_D1_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_PSSK3, SK_"test_setSortedArrQsortrdp_D1_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_PSSK2, SK_"test_setSortedArrQsortrdp_D1_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrQsortrdp_D1_PSSK1, SK_"test_setSortedArrQsortrdp_D1_PSSK1")
#endif
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrBubble_D0_SK5, SK_"test_setSortedArrBubble_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrBubble_D0_SK4, SK_"test_setSortedArrBubble_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrBubble_D0_SK3, SK_"test_setSortedArrBubble_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrBubble_D0_SK2, SK_"test_setSortedArrBubble_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrBubble_D0_SK1, SK_"test_setSortedArrBubble_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrBubble_D1_SK5, SK_"test_setSortedArrBubble_D1_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrBubble_D1_SK4, SK_"test_setSortedArrBubble_D1_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrBubble_D1_SK3, SK_"test_setSortedArrBubble_D1_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrBubble_D1_SK2, SK_"test_setSortedArrBubble_D1_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrBubble_D1_SK1, SK_"test_setSortedArrBubble_D1_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setSortedArrBubble_D1_IK5, SK_"test_setSortedArrBubble_D1_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_setSortedArrBubble_D1_IK4, SK_"test_setSortedArrBubble_D1_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_setSortedArrBubble_D1_IK3, SK_"test_setSortedArrBubble_D1_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_setSortedArrBubble_D1_IK2, SK_"test_setSortedArrBubble_D1_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_setSortedArrBubble_D1_IK1, SK_"test_setSortedArrBubble_D1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setSortedArrBubble_D1_LK5, SK_"test_setSortedArrBubble_D1_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_setSortedArrBubble_D1_LK4, SK_"test_setSortedArrBubble_D1_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_setSortedArrBubble_D1_LK3, SK_"test_setSortedArrBubble_D1_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_setSortedArrBubble_D1_LK2, SK_"test_setSortedArrBubble_D1_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_setSortedArrBubble_D1_LK1, SK_"test_setSortedArrBubble_D1_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setSortedArrBubble_D1_CK5, SK_"test_setSortedArrBubble_D1_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setSortedArrBubble_D1_CK4, SK_"test_setSortedArrBubble_D1_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setSortedArrBubble_D1_CK3, SK_"test_setSortedArrBubble_D1_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setSortedArrBubble_D1_CK2, SK_"test_setSortedArrBubble_D1_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setSortedArrBubble_D1_CK1, SK_"test_setSortedArrBubble_D1_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setSortedArrBubble_D1_RK5, SK_"test_setSortedArrBubble_D1_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setSortedArrBubble_D1_RK4, SK_"test_setSortedArrBubble_D1_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setSortedArrBubble_D1_RK3, SK_"test_setSortedArrBubble_D1_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setSortedArrBubble_D1_RK2, SK_"test_setSortedArrBubble_D1_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setSortedArrBubble_D1_RK1, SK_"test_setSortedArrBubble_D1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
#if     SK5_ENABLED
        call test%run(test_setSortedArrBubble_D1_PSSK5, SK_"test_setSortedArrBubble_D1_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrBubble_D1_PSSK4, SK_"test_setSortedArrBubble_D1_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrBubble_D1_PSSK3, SK_"test_setSortedArrBubble_D1_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrBubble_D1_PSSK2, SK_"test_setSortedArrBubble_D1_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrBubble_D1_PSSK1, SK_"test_setSortedArrBubble_D1_PSSK1")
#endif
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrHeapi_D0_SK5, SK_"test_setSortedArrHeapi_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrHeapi_D0_SK4, SK_"test_setSortedArrHeapi_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrHeapi_D0_SK3, SK_"test_setSortedArrHeapi_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrHeapi_D0_SK2, SK_"test_setSortedArrHeapi_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrHeapi_D0_SK1, SK_"test_setSortedArrHeapi_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrHeapi_D1_SK5, SK_"test_setSortedArrHeapi_D1_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrHeapi_D1_SK4, SK_"test_setSortedArrHeapi_D1_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrHeapi_D1_SK3, SK_"test_setSortedArrHeapi_D1_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrHeapi_D1_SK2, SK_"test_setSortedArrHeapi_D1_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrHeapi_D1_SK1, SK_"test_setSortedArrHeapi_D1_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setSortedArrHeapi_D1_IK5, SK_"test_setSortedArrHeapi_D1_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_setSortedArrHeapi_D1_IK4, SK_"test_setSortedArrHeapi_D1_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_setSortedArrHeapi_D1_IK3, SK_"test_setSortedArrHeapi_D1_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_setSortedArrHeapi_D1_IK2, SK_"test_setSortedArrHeapi_D1_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_setSortedArrHeapi_D1_IK1, SK_"test_setSortedArrHeapi_D1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setSortedArrHeapi_D1_LK5, SK_"test_setSortedArrHeapi_D1_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_setSortedArrHeapi_D1_LK4, SK_"test_setSortedArrHeapi_D1_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_setSortedArrHeapi_D1_LK3, SK_"test_setSortedArrHeapi_D1_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_setSortedArrHeapi_D1_LK2, SK_"test_setSortedArrHeapi_D1_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_setSortedArrHeapi_D1_LK1, SK_"test_setSortedArrHeapi_D1_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setSortedArrHeapi_D1_CK5, SK_"test_setSortedArrHeapi_D1_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setSortedArrHeapi_D1_CK4, SK_"test_setSortedArrHeapi_D1_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setSortedArrHeapi_D1_CK3, SK_"test_setSortedArrHeapi_D1_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setSortedArrHeapi_D1_CK2, SK_"test_setSortedArrHeapi_D1_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setSortedArrHeapi_D1_CK1, SK_"test_setSortedArrHeapi_D1_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setSortedArrHeapi_D1_RK5, SK_"test_setSortedArrHeapi_D1_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setSortedArrHeapi_D1_RK4, SK_"test_setSortedArrHeapi_D1_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setSortedArrHeapi_D1_RK3, SK_"test_setSortedArrHeapi_D1_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setSortedArrHeapi_D1_RK2, SK_"test_setSortedArrHeapi_D1_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setSortedArrHeapi_D1_RK1, SK_"test_setSortedArrHeapi_D1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
#if     SK5_ENABLED
        call test%run(test_setSortedArrHeapi_D1_PSSK5, SK_"test_setSortedArrHeapi_D1_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrHeapi_D1_PSSK4, SK_"test_setSortedArrHeapi_D1_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrHeapi_D1_PSSK3, SK_"test_setSortedArrHeapi_D1_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrHeapi_D1_PSSK2, SK_"test_setSortedArrHeapi_D1_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrHeapi_D1_PSSK1, SK_"test_setSortedArrHeapi_D1_PSSK1")
#endif
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrHeapr_D0_SK5, SK_"test_setSortedArrHeapr_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrHeapr_D0_SK4, SK_"test_setSortedArrHeapr_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrHeapr_D0_SK3, SK_"test_setSortedArrHeapr_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrHeapr_D0_SK2, SK_"test_setSortedArrHeapr_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrHeapr_D0_SK1, SK_"test_setSortedArrHeapr_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrHeapr_D1_SK5, SK_"test_setSortedArrHeapr_D1_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrHeapr_D1_SK4, SK_"test_setSortedArrHeapr_D1_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrHeapr_D1_SK3, SK_"test_setSortedArrHeapr_D1_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrHeapr_D1_SK2, SK_"test_setSortedArrHeapr_D1_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrHeapr_D1_SK1, SK_"test_setSortedArrHeapr_D1_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setSortedArrHeapr_D1_IK5, SK_"test_setSortedArrHeapr_D1_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_setSortedArrHeapr_D1_IK4, SK_"test_setSortedArrHeapr_D1_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_setSortedArrHeapr_D1_IK3, SK_"test_setSortedArrHeapr_D1_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_setSortedArrHeapr_D1_IK2, SK_"test_setSortedArrHeapr_D1_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_setSortedArrHeapr_D1_IK1, SK_"test_setSortedArrHeapr_D1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setSortedArrHeapr_D1_LK5, SK_"test_setSortedArrHeapr_D1_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_setSortedArrHeapr_D1_LK4, SK_"test_setSortedArrHeapr_D1_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_setSortedArrHeapr_D1_LK3, SK_"test_setSortedArrHeapr_D1_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_setSortedArrHeapr_D1_LK2, SK_"test_setSortedArrHeapr_D1_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_setSortedArrHeapr_D1_LK1, SK_"test_setSortedArrHeapr_D1_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setSortedArrHeapr_D1_CK5, SK_"test_setSortedArrHeapr_D1_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setSortedArrHeapr_D1_CK4, SK_"test_setSortedArrHeapr_D1_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setSortedArrHeapr_D1_CK3, SK_"test_setSortedArrHeapr_D1_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setSortedArrHeapr_D1_CK2, SK_"test_setSortedArrHeapr_D1_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setSortedArrHeapr_D1_CK1, SK_"test_setSortedArrHeapr_D1_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setSortedArrHeapr_D1_RK5, SK_"test_setSortedArrHeapr_D1_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setSortedArrHeapr_D1_RK4, SK_"test_setSortedArrHeapr_D1_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setSortedArrHeapr_D1_RK3, SK_"test_setSortedArrHeapr_D1_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setSortedArrHeapr_D1_RK2, SK_"test_setSortedArrHeapr_D1_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setSortedArrHeapr_D1_RK1, SK_"test_setSortedArrHeapr_D1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
#if     SK5_ENABLED
        call test%run(test_setSortedArrHeapr_D1_PSSK5, SK_"test_setSortedArrHeapr_D1_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrHeapr_D1_PSSK4, SK_"test_setSortedArrHeapr_D1_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrHeapr_D1_PSSK3, SK_"test_setSortedArrHeapr_D1_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrHeapr_D1_PSSK2, SK_"test_setSortedArrHeapr_D1_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrHeapr_D1_PSSK1, SK_"test_setSortedArrHeapr_D1_PSSK1")
#endif
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrInsertionl_D0_SK5, SK_"test_setSortedArrInsertionl_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrInsertionl_D0_SK4, SK_"test_setSortedArrInsertionl_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrInsertionl_D0_SK3, SK_"test_setSortedArrInsertionl_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrInsertionl_D0_SK2, SK_"test_setSortedArrInsertionl_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrInsertionl_D0_SK1, SK_"test_setSortedArrInsertionl_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_SK5, SK_"test_setSortedArrInsertionl_D1_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_SK4, SK_"test_setSortedArrInsertionl_D1_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_SK3, SK_"test_setSortedArrInsertionl_D1_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_SK2, SK_"test_setSortedArrInsertionl_D1_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_SK1, SK_"test_setSortedArrInsertionl_D1_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_IK5, SK_"test_setSortedArrInsertionl_D1_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_IK4, SK_"test_setSortedArrInsertionl_D1_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_IK3, SK_"test_setSortedArrInsertionl_D1_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_IK2, SK_"test_setSortedArrInsertionl_D1_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_IK1, SK_"test_setSortedArrInsertionl_D1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_LK5, SK_"test_setSortedArrInsertionl_D1_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_LK4, SK_"test_setSortedArrInsertionl_D1_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_LK3, SK_"test_setSortedArrInsertionl_D1_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_LK2, SK_"test_setSortedArrInsertionl_D1_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_LK1, SK_"test_setSortedArrInsertionl_D1_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_CK5, SK_"test_setSortedArrInsertionl_D1_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_CK4, SK_"test_setSortedArrInsertionl_D1_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_CK3, SK_"test_setSortedArrInsertionl_D1_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_CK2, SK_"test_setSortedArrInsertionl_D1_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_CK1, SK_"test_setSortedArrInsertionl_D1_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_RK5, SK_"test_setSortedArrInsertionl_D1_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_RK4, SK_"test_setSortedArrInsertionl_D1_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_RK3, SK_"test_setSortedArrInsertionl_D1_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_RK2, SK_"test_setSortedArrInsertionl_D1_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_RK1, SK_"test_setSortedArrInsertionl_D1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
#if     SK5_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_PSSK5, SK_"test_setSortedArrInsertionl_D1_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_PSSK4, SK_"test_setSortedArrInsertionl_D1_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_PSSK3, SK_"test_setSortedArrInsertionl_D1_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_PSSK2, SK_"test_setSortedArrInsertionl_D1_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrInsertionl_D1_PSSK1, SK_"test_setSortedArrInsertionl_D1_PSSK1")
#endif
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrInsertionb_D0_SK5, SK_"test_setSortedArrInsertionb_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrInsertionb_D0_SK4, SK_"test_setSortedArrInsertionb_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrInsertionb_D0_SK3, SK_"test_setSortedArrInsertionb_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrInsertionb_D0_SK2, SK_"test_setSortedArrInsertionb_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrInsertionb_D0_SK1, SK_"test_setSortedArrInsertionb_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_SK5, SK_"test_setSortedArrInsertionb_D1_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_SK4, SK_"test_setSortedArrInsertionb_D1_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_SK3, SK_"test_setSortedArrInsertionb_D1_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_SK2, SK_"test_setSortedArrInsertionb_D1_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_SK1, SK_"test_setSortedArrInsertionb_D1_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_IK5, SK_"test_setSortedArrInsertionb_D1_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_IK4, SK_"test_setSortedArrInsertionb_D1_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_IK3, SK_"test_setSortedArrInsertionb_D1_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_IK2, SK_"test_setSortedArrInsertionb_D1_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_IK1, SK_"test_setSortedArrInsertionb_D1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_LK5, SK_"test_setSortedArrInsertionb_D1_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_LK4, SK_"test_setSortedArrInsertionb_D1_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_LK3, SK_"test_setSortedArrInsertionb_D1_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_LK2, SK_"test_setSortedArrInsertionb_D1_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_LK1, SK_"test_setSortedArrInsertionb_D1_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_CK5, SK_"test_setSortedArrInsertionb_D1_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_CK4, SK_"test_setSortedArrInsertionb_D1_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_CK3, SK_"test_setSortedArrInsertionb_D1_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_CK2, SK_"test_setSortedArrInsertionb_D1_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_CK1, SK_"test_setSortedArrInsertionb_D1_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_RK5, SK_"test_setSortedArrInsertionb_D1_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_RK4, SK_"test_setSortedArrInsertionb_D1_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_RK3, SK_"test_setSortedArrInsertionb_D1_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_RK2, SK_"test_setSortedArrInsertionb_D1_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_RK1, SK_"test_setSortedArrInsertionb_D1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
#if     SK5_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_PSSK5, SK_"test_setSortedArrInsertionb_D1_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_PSSK4, SK_"test_setSortedArrInsertionb_D1_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_PSSK3, SK_"test_setSortedArrInsertionb_D1_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_PSSK2, SK_"test_setSortedArrInsertionb_D1_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrInsertionb_D1_PSSK1, SK_"test_setSortedArrInsertionb_D1_PSSK1")
#endif
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrMerger_D0_SK5, SK_"test_setSortedArrMerger_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrMerger_D0_SK4, SK_"test_setSortedArrMerger_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrMerger_D0_SK3, SK_"test_setSortedArrMerger_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrMerger_D0_SK2, SK_"test_setSortedArrMerger_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrMerger_D0_SK1, SK_"test_setSortedArrMerger_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrMerger_D1_SK5, SK_"test_setSortedArrMerger_D1_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrMerger_D1_SK4, SK_"test_setSortedArrMerger_D1_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrMerger_D1_SK3, SK_"test_setSortedArrMerger_D1_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrMerger_D1_SK2, SK_"test_setSortedArrMerger_D1_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrMerger_D1_SK1, SK_"test_setSortedArrMerger_D1_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setSortedArrMerger_D1_IK5, SK_"test_setSortedArrMerger_D1_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_setSortedArrMerger_D1_IK4, SK_"test_setSortedArrMerger_D1_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_setSortedArrMerger_D1_IK3, SK_"test_setSortedArrMerger_D1_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_setSortedArrMerger_D1_IK2, SK_"test_setSortedArrMerger_D1_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_setSortedArrMerger_D1_IK1, SK_"test_setSortedArrMerger_D1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setSortedArrMerger_D1_LK5, SK_"test_setSortedArrMerger_D1_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_setSortedArrMerger_D1_LK4, SK_"test_setSortedArrMerger_D1_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_setSortedArrMerger_D1_LK3, SK_"test_setSortedArrMerger_D1_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_setSortedArrMerger_D1_LK2, SK_"test_setSortedArrMerger_D1_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_setSortedArrMerger_D1_LK1, SK_"test_setSortedArrMerger_D1_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setSortedArrMerger_D1_CK5, SK_"test_setSortedArrMerger_D1_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setSortedArrMerger_D1_CK4, SK_"test_setSortedArrMerger_D1_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setSortedArrMerger_D1_CK3, SK_"test_setSortedArrMerger_D1_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setSortedArrMerger_D1_CK2, SK_"test_setSortedArrMerger_D1_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setSortedArrMerger_D1_CK1, SK_"test_setSortedArrMerger_D1_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setSortedArrMerger_D1_RK5, SK_"test_setSortedArrMerger_D1_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setSortedArrMerger_D1_RK4, SK_"test_setSortedArrMerger_D1_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setSortedArrMerger_D1_RK3, SK_"test_setSortedArrMerger_D1_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setSortedArrMerger_D1_RK2, SK_"test_setSortedArrMerger_D1_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setSortedArrMerger_D1_RK1, SK_"test_setSortedArrMerger_D1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
#if     SK5_ENABLED
        call test%run(test_setSortedArrMerger_D1_PSSK5, SK_"test_setSortedArrMerger_D1_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrMerger_D1_PSSK4, SK_"test_setSortedArrMerger_D1_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrMerger_D1_PSSK3, SK_"test_setSortedArrMerger_D1_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrMerger_D1_PSSK2, SK_"test_setSortedArrMerger_D1_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrMerger_D1_PSSK1, SK_"test_setSortedArrMerger_D1_PSSK1")
#endif
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrSelection_D0_SK5, SK_"test_setSortedArrSelection_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrSelection_D0_SK4, SK_"test_setSortedArrSelection_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrSelection_D0_SK3, SK_"test_setSortedArrSelection_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrSelection_D0_SK2, SK_"test_setSortedArrSelection_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrSelection_D0_SK1, SK_"test_setSortedArrSelection_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrSelection_D1_SK5, SK_"test_setSortedArrSelection_D1_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrSelection_D1_SK4, SK_"test_setSortedArrSelection_D1_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrSelection_D1_SK3, SK_"test_setSortedArrSelection_D1_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrSelection_D1_SK2, SK_"test_setSortedArrSelection_D1_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrSelection_D1_SK1, SK_"test_setSortedArrSelection_D1_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setSortedArrSelection_D1_IK5, SK_"test_setSortedArrSelection_D1_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_setSortedArrSelection_D1_IK4, SK_"test_setSortedArrSelection_D1_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_setSortedArrSelection_D1_IK3, SK_"test_setSortedArrSelection_D1_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_setSortedArrSelection_D1_IK2, SK_"test_setSortedArrSelection_D1_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_setSortedArrSelection_D1_IK1, SK_"test_setSortedArrSelection_D1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setSortedArrSelection_D1_LK5, SK_"test_setSortedArrSelection_D1_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_setSortedArrSelection_D1_LK4, SK_"test_setSortedArrSelection_D1_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_setSortedArrSelection_D1_LK3, SK_"test_setSortedArrSelection_D1_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_setSortedArrSelection_D1_LK2, SK_"test_setSortedArrSelection_D1_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_setSortedArrSelection_D1_LK1, SK_"test_setSortedArrSelection_D1_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setSortedArrSelection_D1_CK5, SK_"test_setSortedArrSelection_D1_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setSortedArrSelection_D1_CK4, SK_"test_setSortedArrSelection_D1_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setSortedArrSelection_D1_CK3, SK_"test_setSortedArrSelection_D1_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setSortedArrSelection_D1_CK2, SK_"test_setSortedArrSelection_D1_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setSortedArrSelection_D1_CK1, SK_"test_setSortedArrSelection_D1_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setSortedArrSelection_D1_RK5, SK_"test_setSortedArrSelection_D1_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setSortedArrSelection_D1_RK4, SK_"test_setSortedArrSelection_D1_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setSortedArrSelection_D1_RK3, SK_"test_setSortedArrSelection_D1_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setSortedArrSelection_D1_RK2, SK_"test_setSortedArrSelection_D1_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setSortedArrSelection_D1_RK1, SK_"test_setSortedArrSelection_D1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
#if     SK5_ENABLED
        call test%run(test_setSortedArrSelection_D1_PSSK5, SK_"test_setSortedArrSelection_D1_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrSelection_D1_PSSK4, SK_"test_setSortedArrSelection_D1_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrSelection_D1_PSSK3, SK_"test_setSortedArrSelection_D1_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrSelection_D1_PSSK2, SK_"test_setSortedArrSelection_D1_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrSelection_D1_PSSK1, SK_"test_setSortedArrSelection_D1_PSSK1")
#endif
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrShell_D0_SK5, SK_"test_setSortedArrShell_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrShell_D0_SK4, SK_"test_setSortedArrShell_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrShell_D0_SK3, SK_"test_setSortedArrShell_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrShell_D0_SK2, SK_"test_setSortedArrShell_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrShell_D0_SK1, SK_"test_setSortedArrShell_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setSortedArrShell_D1_SK5, SK_"test_setSortedArrShell_D1_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrShell_D1_SK4, SK_"test_setSortedArrShell_D1_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrShell_D1_SK3, SK_"test_setSortedArrShell_D1_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrShell_D1_SK2, SK_"test_setSortedArrShell_D1_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrShell_D1_SK1, SK_"test_setSortedArrShell_D1_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setSortedArrShell_D1_IK5, SK_"test_setSortedArrShell_D1_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_setSortedArrShell_D1_IK4, SK_"test_setSortedArrShell_D1_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_setSortedArrShell_D1_IK3, SK_"test_setSortedArrShell_D1_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_setSortedArrShell_D1_IK2, SK_"test_setSortedArrShell_D1_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_setSortedArrShell_D1_IK1, SK_"test_setSortedArrShell_D1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setSortedArrShell_D1_LK5, SK_"test_setSortedArrShell_D1_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_setSortedArrShell_D1_LK4, SK_"test_setSortedArrShell_D1_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_setSortedArrShell_D1_LK3, SK_"test_setSortedArrShell_D1_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_setSortedArrShell_D1_LK2, SK_"test_setSortedArrShell_D1_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_setSortedArrShell_D1_LK1, SK_"test_setSortedArrShell_D1_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setSortedArrShell_D1_CK5, SK_"test_setSortedArrShell_D1_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setSortedArrShell_D1_CK4, SK_"test_setSortedArrShell_D1_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setSortedArrShell_D1_CK3, SK_"test_setSortedArrShell_D1_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setSortedArrShell_D1_CK2, SK_"test_setSortedArrShell_D1_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setSortedArrShell_D1_CK1, SK_"test_setSortedArrShell_D1_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setSortedArrShell_D1_RK5, SK_"test_setSortedArrShell_D1_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setSortedArrShell_D1_RK4, SK_"test_setSortedArrShell_D1_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setSortedArrShell_D1_RK3, SK_"test_setSortedArrShell_D1_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setSortedArrShell_D1_RK2, SK_"test_setSortedArrShell_D1_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setSortedArrShell_D1_RK1, SK_"test_setSortedArrShell_D1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
#if     SK5_ENABLED
        call test%run(test_setSortedArrShell_D1_PSSK5, SK_"test_setSortedArrShell_D1_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setSortedArrShell_D1_PSSK4, SK_"test_setSortedArrShell_D1_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setSortedArrShell_D1_PSSK3, SK_"test_setSortedArrShell_D1_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setSortedArrShell_D1_PSSK2, SK_"test_setSortedArrShell_D1_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setSortedArrShell_D1_PSSK1, SK_"test_setSortedArrShell_D1_PSSK1")
#endif
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !error stop
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if 0
    !>  \brief
    !> Test the performance of different sorting algorithms in [pm_arraySort](@ref pm_arraySort).
    ! LCOV_EXCL_START
    function test_performance_1() result(assertion)

        implicit none
        logical(LK) :: assertion
        character(:, SK), allocatable :: methodName

        assertion = .true._LK
        test%File = test%openFile()

        !write(test%File%unit,"(*(g0,:,' '))") "Smooth sort" ! Smooth search is too slow to request arrays sizes > 10000.
        !call bench(sortSmooth,10000_IK)
        !write(test%File%unit,"(*(g0,:,' '))")

        methodName = "Bubble sort"
        write(test%File%unit,"(*(g0,:,' '))") methodName
        call bench(setSortedBubble_RK2,10000_IK)
        write(test%File%unit,"(*(g0,:,' '))")

        methodName = "Selection sort"
        write(test%File%unit,"(*(g0,:,' '))") methodName
        call bench(setSortedSelection_RK2,10000_IK)
        write(test%File%unit,"(*(g0,:,' '))")

        methodName = "insertion sort"
        write(test%File%unit,"(*(g0,:,' '))") methodName
        call bench(setSortedInsertion_RK2,10000_IK)
        write(test%File%unit,"(*(g0,:,' '))")

        methodName = "Binary insertion sort"
        write(test%File%unit,"(*(g0,:,' '))") methodName
        call bench(setSortedInsertionb_RK2,10000_IK)
        write(test%File%unit,"(*(g0,:,' '))")

        !methodName = "Pair insertion sort"
        !write(test%File%unit,"(*(g0,:,' '))") methodName
        !call bench(setSortedInsertionPair,10000_IK)
        !write(test%File%unit,"(*(g0,:,' '))")

        methodName = "Double insertion sort"
        write(test%File%unit,"(*(g0,:,' '))") methodName
        call bench(setSortedInsertionDouble_RK2,10000_IK)
        write(test%File%unit,"(*(g0,:,' '))")

        methodName = "Shell sort"
        write(test%File%unit,"(*(g0,:,' '))") methodName
        call bench(setSortedShell_RK2,100000_IK)
        write(test%File%unit,"(*(g0,:,' '))")

        methodName = "Merge sort"
        write(test%File%unit,"(*(g0,:,' '))") methodName
        call bench(setSortedMergeRecDefCom_RK2,100000_IK)
        write(test%File%unit,"(*(g0,:,' '))")

        methodName = "Heap sort"
        write(test%File%unit,"(*(g0,:,' '))") methodName
        call bench(setSortedHeap_RK2,100000_IK)
        write(test%File%unit,"(*(g0,:,' '))")

        methodName = "Heap sort recursive"
        write(test%File%unit,"(*(g0,:,' '))") methodName
        call bench(setSortedHeapRec_RK2,100000_IK)
        write(test%File%unit,"(*(g0,:,' '))")

        methodName = "Quicksort"
        write(test%File%unit,"(*(g0,:,' '))") methodName
        call bench(setSortedRec_RK2,100000_IK)
        write(test%File%unit,"(*(g0,:,' '))")

        methodName = "Dual pivot quicksort"
        write(test%File%unit,"(*(g0,:,' '))") methodName
        call bench(setSortedDualPivot_RK2,100000_IK)
        write(test%File%unit,"(*(g0,:,' '))")

        methodName = "Quicksort (no recursion)"
        write(test%File%unit,"(*(g0,:,' '))") methodName
        call bench(sortWrapper,100000_IK)
        write(test%File%unit,"(*(g0,:,' '))")

        close(test%File%unit)

    contains

        subroutine sortWrapper(Array)
            real(RK), intent(inout), contiguous :: Array(:)
            call setSorted(Array)
        end subroutine

        subroutine bench(sort,arraySize)
            abstract interface
                subroutine setSorted_proc_RK2(Array)
                    use pm_kind, only: RK2
                    real(RK2)  , intent(inout), contiguous :: Array(:)
                end subroutine setSorted_proc_RK2
            end interface
            procedure(setSorted_proc_RK2)   :: sort
            real(RK)    , allocatable   :: Array(:)
            integer(IK) , intent(in)    :: arraySize
            integer(IK)                 :: i,j,off,t1,t2,rate
            real(RK)                    :: check
            allocate(Array(arraySize))
            i = 10_IK
            ! find rate
            call system_clock(t1,rate)
            write(test%File%unit,"(*(g0,:,' '))") "      Size      time/element"
            do while(i <= arraySize)
                call random_number(Array)
                check = sum(Array(1:i))
                off = 1
                call system_clock(t1)
                do while(off+i-1 <= arraySize)
                    call setSorted(Array(off:off+i-1))
                    off = off + i
                end do
                call system_clock(t2)
                ! Check result is ordered
                do j = 2, i
                    if (Array(j) < Array(j-1)) then
                        assertion = .false._LK
                        write(test%File%unit,"(*(g0,:,' '))") "Sort error at ",j," value ",Array(j)
                        if (i <= 100_IK) then
                            write(test%File%unit,"(*(g0,:,' '))") Array(1:i)
                        else
                            !write(test%File%unit,"(*(g0,:,' '))") Array(j-3:j+3)
                            write(test%File%unit,"(*(g0,:,' '))") Array(j-1), Array(j)
                        end if
                        write(*,*) methodName
                        !return
                    end if
                end do
                ! Check sum of result is approximately the same as it was
                check = check - sum(Array(1:i))
                if (check > 1e-15_RK * i*sqrt(real(i,RK))) then
                    write(test%File%unit,"(*(g0,:,' '))") "Large difference in sums before and after",check
                end if

                write(test%File%unit,"(*(g0,:,' '))")i,real(t2-t1,RK)/(real(rate,RK)*i*(arraySize/i))
                i = 10 * i
            end do
            deallocate(Array)
        end subroutine bench

    end function test_performance_1
    ! LCOV_EXCL_STOP
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_arraySort