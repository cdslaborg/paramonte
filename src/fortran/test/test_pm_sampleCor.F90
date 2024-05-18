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
!>  This module contains tests of the module [pm_sampleCor](@ref pm_sampleCor).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Thursday 12:45 AM, August 20, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_sampleCor

    use pm_sampleCor
    use pm_kind, only: IK
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! getCor

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_getCor_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_getCor_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_getCor_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_getCor_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_getCor_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getCor_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getCor_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getCor_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getCor_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getCor_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! setCor

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setCor_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setCor_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setCor_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setCor_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setCor_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setCor_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setCor_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setCor_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setCor_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setCor_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! getRho

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_getRho_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_getRho_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_getRho_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_getRho_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_getRho_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_getRho_XY_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_getRho_XY_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_getRho_XY_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_getRho_XY_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_getRho_XY_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_getRho_XY_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_getRho_XY_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_getRho_XY_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_getRho_XY_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_getRho_XY_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getRho_XY_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getRho_XY_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getRho_XY_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getRho_XY_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getRho_XY_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_getRho_XY_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_getRho_XY_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_getRho_XY_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_getRho_XY_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_getRho_XY_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        module function test_getRho_XY_BSSK() result(assertion); logical(LK) :: assertion; end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_getRho_D2_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_getRho_D2_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_getRho_D2_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_getRho_D2_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_getRho_D2_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_getRho_D2_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_getRho_D2_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_getRho_D2_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_getRho_D2_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_getRho_D2_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getRho_D2_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getRho_D2_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getRho_D2_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getRho_D2_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getRho_D2_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_getRho_D2_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_getRho_D2_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_getRho_D2_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_getRho_D2_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_getRho_D2_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        module function test_getRho_D2_BSSK() result(assertion); logical(LK) :: assertion; end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! setRho

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setRho_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setRho_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setRho_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setRho_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setRho_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setRho_XY_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setRho_XY_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setRho_XY_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setRho_XY_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setRho_XY_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setRho_XY_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setRho_XY_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setRho_XY_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setRho_XY_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setRho_XY_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setRho_XY_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setRho_XY_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setRho_XY_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setRho_XY_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setRho_XY_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setRho_XY_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setRho_XY_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setRho_XY_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setRho_XY_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setRho_XY_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        module function test_setRho_XY_BSSK() result(assertion); logical(LK) :: assertion; end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setRho_D2_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setRho_D2_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setRho_D2_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setRho_D2_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setRho_D2_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setRho_D2_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setRho_D2_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setRho_D2_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setRho_D2_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setRho_D2_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setRho_D2_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setRho_D2_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setRho_D2_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setRho_D2_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setRho_D2_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setRho_D2_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setRho_D2_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setRho_D2_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setRho_D2_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setRho_D2_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        module function test_setRho_D2_BSSK() result(assertion); logical(LK) :: assertion; end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! setCordance

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setCordance_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setCordance_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setCordance_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setCordance_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setCordance_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setCordance_D1_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setCordance_D1_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setCordance_D1_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setCordance_D1_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setCordance_D1_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setCordance_D1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setCordance_D1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setCordance_D1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setCordance_D1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setCordance_D1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setCordance_D1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setCordance_D1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setCordance_D1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setCordance_D1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setCordance_D1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setCordance_D1_PSSK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setCordance_D1_PSSK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setCordance_D1_PSSK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setCordance_D1_PSSK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setCordance_D1_PSSK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        module function test_setCordance_D1_BSSK() result(assertion); logical(LK) :: assertion; end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none

        test = test_type(MODULE_NAME)

        ! getCor

        block

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_getCor_CK5, SK_"test_getCor_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_getCor_CK4, SK_"test_getCor_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_getCor_CK3, SK_"test_getCor_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_getCor_CK2, SK_"test_getCor_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_getCor_CK1, SK_"test_getCor_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getCor_RK5, SK_"test_getCor_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getCor_RK4, SK_"test_getCor_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getCor_RK3, SK_"test_getCor_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getCor_RK2, SK_"test_getCor_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getCor_RK1, SK_"test_getCor_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end block

        ! setCor

        block

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setCor_CK5, SK_"test_setCor_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setCor_CK4, SK_"test_setCor_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setCor_CK3, SK_"test_setCor_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setCor_CK2, SK_"test_setCor_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setCor_CK1, SK_"test_setCor_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setCor_RK5, SK_"test_setCor_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setCor_RK4, SK_"test_setCor_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setCor_RK3, SK_"test_setCor_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setCor_RK2, SK_"test_setCor_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setCor_RK1, SK_"test_setCor_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end block

        ! getRho

        block

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_getRho_D0_SK5, SK_"test_getRho_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_getRho_D0_SK4, SK_"test_getRho_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_getRho_D0_SK3, SK_"test_getRho_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_getRho_D0_SK2, SK_"test_getRho_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_getRho_D0_SK1, SK_"test_getRho_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_getRho_XY_SK5, SK_"test_getRho_XY_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_getRho_XY_SK4, SK_"test_getRho_XY_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_getRho_XY_SK3, SK_"test_getRho_XY_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_getRho_XY_SK2, SK_"test_getRho_XY_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_getRho_XY_SK1, SK_"test_getRho_XY_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_getRho_XY_IK5, SK_"test_getRho_XY_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_getRho_XY_IK4, SK_"test_getRho_XY_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_getRho_XY_IK3, SK_"test_getRho_XY_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_getRho_XY_IK2, SK_"test_getRho_XY_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_getRho_XY_IK1, SK_"test_getRho_XY_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getRho_XY_RK5, SK_"test_getRho_XY_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getRho_XY_RK4, SK_"test_getRho_XY_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getRho_XY_RK3, SK_"test_getRho_XY_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getRho_XY_RK2, SK_"test_getRho_XY_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getRho_XY_RK1, SK_"test_getRho_XY_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     PDT_ENABLED
#if     SK5_ENABLED
        call test%run(test_getRho_XY_PSSK5, SK_"test_getRho_XY_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_getRho_XY_PSSK4, SK_"test_getRho_XY_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_getRho_XY_PSSK3, SK_"test_getRho_XY_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_getRho_XY_PSSK2, SK_"test_getRho_XY_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_getRho_XY_PSSK1, SK_"test_getRho_XY_PSSK1")
#endif
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%run(test_getRho_XY_BSSK, SK_"test_getRho_XY_BSSK")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_getRho_D2_SK5, SK_"test_getRho_D2_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_getRho_D2_SK4, SK_"test_getRho_D2_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_getRho_D2_SK3, SK_"test_getRho_D2_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_getRho_D2_SK2, SK_"test_getRho_D2_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_getRho_D2_SK1, SK_"test_getRho_D2_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_getRho_D2_IK5, SK_"test_getRho_D2_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_getRho_D2_IK4, SK_"test_getRho_D2_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_getRho_D2_IK3, SK_"test_getRho_D2_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_getRho_D2_IK2, SK_"test_getRho_D2_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_getRho_D2_IK1, SK_"test_getRho_D2_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getRho_D2_RK5, SK_"test_getRho_D2_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getRho_D2_RK4, SK_"test_getRho_D2_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getRho_D2_RK3, SK_"test_getRho_D2_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getRho_D2_RK2, SK_"test_getRho_D2_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getRho_D2_RK1, SK_"test_getRho_D2_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     PDT_ENABLED
#if     SK5_ENABLED
        call test%run(test_getRho_D2_PSSK5, SK_"test_getRho_D2_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_getRho_D2_PSSK4, SK_"test_getRho_D2_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_getRho_D2_PSSK3, SK_"test_getRho_D2_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_getRho_D2_PSSK2, SK_"test_getRho_D2_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_getRho_D2_PSSK1, SK_"test_getRho_D2_PSSK1")
#endif
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%run(test_getRho_D2_BSSK, SK_"test_getRho_D2_BSSK")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end block

        ! setRho

        block

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setRho_D0_SK5, SK_"test_setRho_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setRho_D0_SK4, SK_"test_setRho_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setRho_D0_SK3, SK_"test_setRho_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setRho_D0_SK2, SK_"test_setRho_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setRho_D0_SK1, SK_"test_setRho_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setRho_XY_SK5, SK_"test_setRho_XY_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setRho_XY_SK4, SK_"test_setRho_XY_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setRho_XY_SK3, SK_"test_setRho_XY_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setRho_XY_SK2, SK_"test_setRho_XY_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setRho_XY_SK1, SK_"test_setRho_XY_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setRho_XY_IK5, SK_"test_setRho_XY_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_setRho_XY_IK4, SK_"test_setRho_XY_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_setRho_XY_IK3, SK_"test_setRho_XY_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_setRho_XY_IK2, SK_"test_setRho_XY_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_setRho_XY_IK1, SK_"test_setRho_XY_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setRho_XY_RK5, SK_"test_setRho_XY_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setRho_XY_RK4, SK_"test_setRho_XY_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setRho_XY_RK3, SK_"test_setRho_XY_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setRho_XY_RK2, SK_"test_setRho_XY_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setRho_XY_RK1, SK_"test_setRho_XY_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     PDT_ENABLED
#if     SK5_ENABLED
        call test%run(test_setRho_XY_PSSK5, SK_"test_setRho_XY_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setRho_XY_PSSK4, SK_"test_setRho_XY_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setRho_XY_PSSK3, SK_"test_setRho_XY_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setRho_XY_PSSK2, SK_"test_setRho_XY_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setRho_XY_PSSK1, SK_"test_setRho_XY_PSSK1")
#endif
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%run(test_setRho_XY_BSSK, SK_"test_setRho_XY_BSSK")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setRho_D2_SK5, SK_"test_setRho_D2_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setRho_D2_SK4, SK_"test_setRho_D2_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setRho_D2_SK3, SK_"test_setRho_D2_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setRho_D2_SK2, SK_"test_setRho_D2_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setRho_D2_SK1, SK_"test_setRho_D2_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setRho_D2_IK5, SK_"test_setRho_D2_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_setRho_D2_IK4, SK_"test_setRho_D2_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_setRho_D2_IK3, SK_"test_setRho_D2_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_setRho_D2_IK2, SK_"test_setRho_D2_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_setRho_D2_IK1, SK_"test_setRho_D2_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setRho_D2_RK5, SK_"test_setRho_D2_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setRho_D2_RK4, SK_"test_setRho_D2_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setRho_D2_RK3, SK_"test_setRho_D2_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setRho_D2_RK2, SK_"test_setRho_D2_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setRho_D2_RK1, SK_"test_setRho_D2_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     PDT_ENABLED
#if     SK5_ENABLED
        call test%run(test_setRho_D2_PSSK5, SK_"test_setRho_D2_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setRho_D2_PSSK4, SK_"test_setRho_D2_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setRho_D2_PSSK3, SK_"test_setRho_D2_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setRho_D2_PSSK2, SK_"test_setRho_D2_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setRho_D2_PSSK1, SK_"test_setRho_D2_PSSK1")
#endif
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%run(test_setRho_D2_BSSK, SK_"test_setRho_D2_BSSK")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end block

        ! setCordance

        block

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setCordance_D0_SK5, SK_"test_setCordance_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setCordance_D0_SK4, SK_"test_setCordance_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setCordance_D0_SK3, SK_"test_setCordance_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setCordance_D0_SK2, SK_"test_setCordance_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setCordance_D0_SK1, SK_"test_setCordance_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setCordance_D1_SK5, SK_"test_setCordance_D1_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setCordance_D1_SK4, SK_"test_setCordance_D1_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setCordance_D1_SK3, SK_"test_setCordance_D1_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setCordance_D1_SK2, SK_"test_setCordance_D1_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setCordance_D1_SK1, SK_"test_setCordance_D1_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setCordance_D1_IK5, SK_"test_setCordance_D1_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_setCordance_D1_IK4, SK_"test_setCordance_D1_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_setCordance_D1_IK3, SK_"test_setCordance_D1_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_setCordance_D1_IK2, SK_"test_setCordance_D1_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_setCordance_D1_IK1, SK_"test_setCordance_D1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setCordance_D1_RK5, SK_"test_setCordance_D1_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setCordance_D1_RK4, SK_"test_setCordance_D1_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setCordance_D1_RK3, SK_"test_setCordance_D1_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setCordance_D1_RK2, SK_"test_setCordance_D1_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setCordance_D1_RK1, SK_"test_setCordance_D1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     PDT_ENABLED
#if     SK5_ENABLED
        call test%run(test_setCordance_D1_PSSK5, SK_"test_setCordance_D1_PSSK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setCordance_D1_PSSK4, SK_"test_setCordance_D1_PSSK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setCordance_D1_PSSK3, SK_"test_setCordance_D1_PSSK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setCordance_D1_PSSK2, SK_"test_setCordance_D1_PSSK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setCordance_D1_PSSK1, SK_"test_setCordance_D1_PSSK1")
#endif
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%run(test_setCordance_D1_BSSK, SK_"test_setCordance_D1_BSSK")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end block

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_sampleCor