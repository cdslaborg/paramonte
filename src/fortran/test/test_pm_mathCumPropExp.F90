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
!>  This module contains tests of the module [pm_mathCumPropExp](@ref pm_mathCumPropExp).
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 12:01 AM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_mathCumPropExp

    use pm_mathCumPropExp ! LCOV_EXCL_LINE
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    use pm_kind, only: LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_getCumPropExpSel_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getCumPropExpSel_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getCumPropExpSel_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getCumPropExpSel_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getCumPropExpSel_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_getCumPropExpSeq_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getCumPropExpSeq_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getCumPropExpSeq_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getCumPropExpSeq_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getCumPropExpSeq_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_setCumPropExpSel_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setCumPropExpSel_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setCumPropExpSel_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setCumPropExpSel_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setCumPropExpSel_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_setCumPropExpSeq_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setCumPropExpSeq_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setCumPropExpSeq_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setCumPropExpSeq_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setCumPropExpSeq_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getCumPropExpSel_RK5, SK_"test_getCumPropExpSel_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getCumPropExpSel_RK4, SK_"test_getCumPropExpSel_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getCumPropExpSel_RK3, SK_"test_getCumPropExpSel_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getCumPropExpSel_RK2, SK_"test_getCumPropExpSel_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getCumPropExpSel_RK1, SK_"test_getCumPropExpSel_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getCumPropExpSeq_RK5, SK_"test_getCumPropExpSeq_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getCumPropExpSeq_RK4, SK_"test_getCumPropExpSeq_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getCumPropExpSeq_RK3, SK_"test_getCumPropExpSeq_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getCumPropExpSeq_RK2, SK_"test_getCumPropExpSeq_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getCumPropExpSeq_RK1, SK_"test_getCumPropExpSeq_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setCumPropExpSel_RK5, SK_"test_setCumPropExpSel_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setCumPropExpSel_RK4, SK_"test_setCumPropExpSel_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setCumPropExpSel_RK3, SK_"test_setCumPropExpSel_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setCumPropExpSel_RK2, SK_"test_setCumPropExpSel_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setCumPropExpSel_RK1, SK_"test_setCumPropExpSel_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setCumPropExpSeq_RK5, SK_"test_setCumPropExpSeq_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setCumPropExpSeq_RK4, SK_"test_setCumPropExpSeq_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setCumPropExpSeq_RK3, SK_"test_setCumPropExpSeq_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setCumPropExpSeq_RK2, SK_"test_setCumPropExpSeq_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setCumPropExpSeq_RK1, SK_"test_setCumPropExpSeq_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%run(test_setCumPropExp_RK, SK_"test_setCumPropExp_RK")

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Test the special legacy interface: genLogCumProp_RK()
    function test_setCumPropExp_RK() result(assertion)

        use pm_kind, only: IK, RK
        use pm_mathCumSum, only: getCumSum ! LCOV_EXCL_LINE
        use pm_distUnif, only: getUnifRand ! LCOV_EXCL_LINE

        logical(LK)                 :: assertion
        real(RK)    , parameter     :: TOL = epsilon(1._RK) * 10._RK
        real(RK)    , allocatable   :: Array(:)
        real(RK)    , allocatable   :: cumPropExp_ref(:)
        real(RK)    , allocatable   :: cumPropExp(:)
        real(RK)    , allocatable   :: diff(:)
        integer(IK)                 :: i

        do i = 1, 200
            Array = getUnifRand(minexponent(0._RK), maxexponent(0._RK), getUnifRand(1_IK, 10_IK))
            cumPropExp_ref = getCumSum(exp(Array - maxval(Array)))
            cumPropExp_ref = cumPropExp_ref / cumPropExp_ref(size(cumPropExp_ref))
            cumPropExp = getCumPropExp_RK(Array, maxval(Array), size(Array, kind = IK))
            diff = abs(cumPropExp - cumPropExp_ref)
            assertion = all(diff < TOL)
            call test%assert(assertion, desc = "The genLogCumProp_RK() must correctly compute CumSumProp for an input array of real of default kind.")
        end do

    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_mathCumPropExp ! LCOV_EXCL_LINE