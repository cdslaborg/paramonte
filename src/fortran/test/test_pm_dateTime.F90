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
!>  This include file contains procedure implementations of the tests of [pm_dateTime](@ref pm_dateTime).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, March 22, 2012, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_dateTime

    use pm_dateTime
    use pm_kind, only: LK
    use pm_test, only: test_type

    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
        module function test_constructDateTimeInt() result(assertion); logical(LK) :: assertion; end function
        module function test_constructDateTimeStr() result(assertion); logical(LK) :: assertion; end function
        module function test_getDateTimeShifted() result(assertion); logical(LK) :: assertion; end function
        module function test_getDateTimeNewZone() result(assertion); logical(LK) :: assertion; end function
        module function test_getCountLeapYears() result(assertion); logical(LK) :: assertion; end function
        module function test_isLastDayInMonth() result(assertion); logical(LK) :: assertion; end function
        module function test_getDateTimeDiff() result(assertion); logical(LK) :: assertion; end function
        module function test_isValidDateTime() result(assertion); logical(LK) :: assertion; end function
        module function test_getMillisecond() result(assertion); logical(LK) :: assertion; end function
        module function test_getDateTimeUTC() result(assertion); logical(LK) :: assertion; end function
        module function test_getDateTime_1() result(assertion); logical(LK) :: assertion; end function
        module function test_getDateTime_2() result(assertion); logical(LK) :: assertion; end function
        module function test_getDateBefore() result(assertion); logical(LK) :: assertion; end function
        module function test_getOrdinalDay() result(assertion); logical(LK) :: assertion; end function
        module function test_getWeekNumber() result(assertion); logical(LK) :: assertion; end function
        module function test_getCountWeeks() result(assertion); logical(LK) :: assertion; end function
        module function test_getWeekDayISO() result(assertion); logical(LK) :: assertion; end function
        module function test_getDateAfter() result(assertion); logical(LK) :: assertion; end function
        module function test_getCountDays() result(assertion); logical(LK) :: assertion; end function
        module function test_getJulianDay() result(assertion); logical(LK) :: assertion; end function
        module function test_getWeekDate() result(assertion); logical(LK) :: assertion; end function
        module function test_isValidZone() result(assertion); logical(LK) :: assertion; end function
        module function test_getWeekYear() result(assertion); logical(LK) :: assertion; end function
        module function test_getWeekDay() result(assertion); logical(LK) :: assertion; end function
        module function test_isMorning() result(assertion); logical(LK) :: assertion; end function
        module function test_getZone() result(assertion); logical(LK) :: assertion; end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none

        test = test_type(MODULE_NAME)
        call test%run(test_constructDateTimeInt, SK_"test_constructDateTimeInt")
        call test%run(test_constructDateTimeStr, SK_"test_constructDateTimeStr")
        call test%run(test_getDateTimeShifted, SK_"test_getDateTimeShifted")
        call test%run(test_getDateTimeNewZone, SK_"test_getDateTimeNewZone")
        call test%run(test_getCountLeapYears, SK_"test_getCountLeapYears")
        call test%run(test_isLastDayInMonth, SK_"test_isLastDayInMonth")
        call test%run(test_getDateTimeDiff, SK_"test_getDateTimeDiff")
        call test%run(test_isValidDateTime, SK_"test_isValidDateTime")
        call test%run(test_getMillisecond, SK_"test_getMillisecond")
        call test%run(test_getDateTimeUTC, SK_"test_getDateTimeUTC")
        call test%run(test_getDateTime_1, SK_"test_getDateTime_1")
        call test%run(test_getDateTime_2, SK_"test_getDateTime_2")
        call test%run(test_getDateAfter, SK_"test_getDateAfter")
        call test%run(test_getDateBefore, SK_"test_getDateBefore")
        call test%run(test_getWeekNumber, SK_"test_getWeekNumber")
        call test%run(test_getCountWeeks, SK_"test_getCountWeeks")
        call test%run(test_getWeekDayISO, SK_"test_getWeekDayISO")
        call test%run(test_getOrdinalDay, SK_"test_getOrdinalDay")
        call test%run(test_getCountDays, SK_"test_getCountDays")
        call test%run(test_getJulianDay, SK_"test_getJulianDay")
        call test%run(test_getWeekDate, SK_"test_getWeekDate")
        call test%run(test_getWeekYear, SK_"test_getWeekYear")
        call test%run(test_isValidZone, SK_"test_isValidZone")
        call test%run(test_getWeekDay, SK_"test_getWeekDay")
        call test%run(test_isMorning, SK_"test_isMorning")
        call test%run(test_getZone, SK_"test_getZone")
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    function test_queryDateTime() result(assertion)
!
!        implicit none
!        logical(LK)         :: assertion
!        type(DateTime_type) :: DateTime
!
!        assertion = .true._LK
!
!        call DateTime%query()
!
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0,:,' '))")
!            write(test%disp%unit,"(*(g0,:,' '))") "           date:", DateTime%date
!            write(test%disp%unit,"(*(g0,:,' '))") "           time:", DateTime%time
!            write(test%disp%unit,"(*(g0,:,' '))") "           zone:", DateTime%zone
!            write(test%disp%unit,"(*(g0,:,' '))") "        century:", DateTime%century
!            write(test%disp%unit,"(*(g0,:,' '))") "           year:", DateTime%year
!            write(test%disp%unit,"(*(g0,:,' '))") "          month:", DateTime%month
!            write(test%disp%unit,"(*(g0,:,' '))") "            day:", DateTime%day
!            write(test%disp%unit,"(*(g0,:,' '))") "           hour:", DateTime%hour
!            write(test%disp%unit,"(*(g0,:,' '))") "         minute:", DateTime%minute
!            write(test%disp%unit,"(*(g0,:,' '))") "         second:", DateTime%second
!            write(test%disp%unit,"(*(g0,:,' '))") "    millisecond:", DateTime%millisecond
!            write(test%disp%unit,"(*(g0,:,' '))") "         values:", DateTime%vals
!            write(test%disp%unit,"(*(g0,:,' '))") "fancyStyleBasic:", DateTime%fancyStyleBasic
!            write(test%disp%unit,"(*(g0,:,' '))") "     fancyStyle:", DateTime%fancyStyle
!            write(test%disp%unit,"(*(g0,:,' '))") 
!            ! LCOV_EXCL_STOP
!        end if
!
!    end function test_queryDateTime
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    ! Test the equivalence of [getNiceDateTime()](@ref pm_dateTime::getnicedatetime) with
!    ! the output of the `fancyStyleBasic()` method of [DateTime_type](@ref pm_dateTime::dateTime_type).
!    function test_getNiceDateTime() result(assertion)
!
!        implicit none
!        logical(LK)                     :: assertion
!        type(DateTime_type)             :: DateTime
!        character(:, SK), allocatable   :: niceDateTime
!
!        call DateTime%query()
!        niceDateTime = getNiceDateTime()
!        assertion = niceDateTime == DateTime%fancyStyleBasic
!
!        if (test%traceable .and. .not. assertion) then
!        ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0,:,' '))")
!            write(test%disp%unit,"(*(g0,:,' '))") "       getNiceDateTime():", niceDateTime
!            write(test%disp%unit,"(*(g0,:,' '))") "DateTime%fancyStyleBasic:", DateTime%fancyStyleBasic
!            write(test%disp%unit,"(*(g0,:,' '))")
!        end if
!        ! LCOV_EXCL_STOP
!
!#if !CODECOV_ENABLED
!        assertion = .true._LK
!#endif
!    end function test_getNiceDateTime

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_dateTime ! LCOV_EXCL_LINE