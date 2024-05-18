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
!>  This module contains tests of the module [pm_str](@ref pm_str).
!>
!>  \final
!>
!>  \author
!>  Amir Shahmoradi

module test_pm_str

    use pm_str
    use pm_kind, only: LK
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
    module function test_getMaxLoc_SK_1 () result(assertion); logical(LK) :: assertion; end function
    module function test_getMaxVal_SK_1 () result(assertion); logical(LK) :: assertion; end function
    module function test_getMinLoc_SK_1 () result(assertion); logical(LK) :: assertion; end function
    module function test_getMinVal_SK_1 () result(assertion); logical(LK) :: assertion; end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none

        test = test_type(MODULE_NAME)

        call test%run(test_getMaxLoc_SK_1, SK_"test_getMaxLoc_SK_1")
        call test%run(test_getMaxVal_SK_1, SK_"test_getMaxVal_SK_1")
        call test%run(test_getMinLoc_SK_1, SK_"test_getMinLoc_SK_1")
        call test%run(test_getMinVal_SK_1, SK_"test_getMinVal_SK_1")
        call test%run(test_getCharVec_SK_1, SK_"test_getCharVec_SK_1")
        call test%run(test_getString_SK_1, SK_"test_getString_SK_1")

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCharVec_SK_1() result(assertion)
        logical(LK) :: assertion
        assertion = all(getCharVec(SK_"ParaMonte") == [character(1,SK) :: 'P','a','r','a','M','o','n','t','e'])
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getString_SK_1() result(assertion)
        logical(LK) :: assertion
        assertion = getCharSeq([character(1,SK) :: 'P','a','r','a','M','o','n','t','e']) == SK_"ParaMonte"
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    function test_padString_1() result(assertion)
!        use pm_kind, only: IK
!        implicit none
!        logical(LK)                     :: assertion
!        integer(IK) , parameter         :: lenPadded = 30_IK
!        character(*, SK), parameter     :: string_nonPadded = "ParaMonte"
!        character(*, SK), parameter     :: stringPadded_ref = "ParaMonte....................."
!        character(*, SK), parameter     :: fill = "."
!        character(:, SK), allocatable   :: strPadded
!        strPadded = getPadded(string_nonPadded, fill, lenPadded)
!        assertion = strPadded == stringPadded_ref .and. len(strPadded) == len(stringPadded_ref)
!        if (test%traceable .and. .not. assertion) then
!        ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "string_nonPadded  : '", string_nonPadded, "'"
!            write(test%disp%unit,"(*(g0))") "stringPadded_ref  : '", stringPadded_ref, "'"
!            write(test%disp%unit,"(*(g0))") "strPadded         : '", strPadded, "'"
!            write(test%disp%unit,"(*(g0))")
!        end if
!        ! LCOV_EXCL_STOP
!    end function test_padString_1
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_padString_2() result(assertion)
!        use pm_kind, only: IK
!        implicit none
!        logical(LK)                     :: assertion
!        integer(IK) , parameter         :: lenPadded = 9_IK
!        character(*, SK), parameter     :: string_nonPadded = "ParaMonte"
!        character(*, SK), parameter     :: stringPadded_ref = "ParaMonte"
!        character(*, SK), parameter     :: fill = "."
!        character(:, SK), allocatable   :: strPadded
!        strPadded = getPadded(string_nonPadded, fill, lenPadded)
!        assertion = strPadded == stringPadded_ref .and. len(strPadded) == len(stringPadded_ref)
!        if (test%traceable .and. .not. assertion) then
!        ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "string_nonPadded  : '", string_nonPadded, "'"
!            write(test%disp%unit,"(*(g0))") "stringPadded_ref  : '", stringPadded_ref, "'"
!            write(test%disp%unit,"(*(g0))") "strPadded         : '", strPadded, "'"
!            write(test%disp%unit,"(*(g0))")
!        end if
!        ! LCOV_EXCL_STOP
!    end function test_padString_2
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !>  \brief
!    !> When `len(string) > lenPadded`, the full string must be returned without any padding.
!    function test_padString_3() result(assertion)
!        use pm_kind, only: IK
!        implicit none
!        logical(LK)                     :: assertion
!        integer(IK) , parameter         :: lenPadded = 5_IK
!        character(*, SK), parameter     :: string_nonPadded = "ParaMonte"
!        character(*, SK), parameter     :: stringPadded_ref = "ParaM"
!        character(*, SK), parameter     :: fill = "."
!        character(:, SK), allocatable   :: strPadded
!        strPadded = getPadded(string_nonPadded, fill, lenPadded)
!        assertion = strPadded == stringPadded_ref .and. len(strPadded) == lenPadded .and. strPadded == stringPadded_ref
!        if (test%traceable .and. .not. assertion) then
!        ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "string_nonPadded  : '", string_nonPadded, "'"
!            write(test%disp%unit,"(*(g0))") "stringPadded_ref  : '", stringPadded_ref, "'"
!            write(test%disp%unit,"(*(g0))") "strPadded         : '", strPadded, "'"
!            write(test%disp%unit,"(*(g0))")
!        end if
!        ! LCOV_EXCL_STOP
!    end function test_padString_3
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_str ! LCOV_EXCL_LINE