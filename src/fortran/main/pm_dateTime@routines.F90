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
!>  This file contains procedure implementations of [pm_dateTime](@ref pm_dateTime).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, March 23, 2012, 5:33 PM, National Institute for Fusion Studies, The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_dateTime) routines

    use pm_val2str, only: getStr
#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isValidZone_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isValidZone_IK_ENABLED 1

    module procedure isValidZone_IK
        use pm_kind, only: IKG => IK
#include "pm_dateTime@routines.inc.F90"
    end procedure

!#if IK5_ENABLED
!    module procedure isValidZone_IK5
!        use pm_kind, only: IKG => IK5
!#include "pm_dateTime@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure isValidZone_IK4
!        use pm_kind, only: IKG => IK4
!#include "pm_dateTime@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure isValidZone_IK3
!        use pm_kind, only: IKG => IK3
!#include "pm_dateTime@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure isValidZone_IK2
!        use pm_kind, only: IKG => IK2
!#include "pm_dateTime@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure isValidZone_IK1
!        use pm_kind, only: IKG => IK1
!#include "pm_dateTime@routines.inc.F90"
!    end procedure
!#endif

#undef isValidZone_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isValidZone_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getZoneAbbrC
        abbr = getZoneAbbr(getZone())
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getZoneAbbrZ
        use pm_arraySort, only: isAscending ! LCOV_EXCL_LINE
        use pm_arraySearch, only: getBin
        use pm_kind, only: SKG => SK
        integer(IK) :: bin
        CHECK_ASSERTION(__LINE__, isAscending(timeZone%zone), SK_"@getZoneAbbr(): The condition `isAscending(timeZone%zone)` must hold. timeZone%zone = "//getStr(timeZone%zone)) ! fpp
        bin = getBin(timeZone%zone, zone) ! we cannot use `findloc()` intrinsic because the time and zone may not be in the constant timezone.
        if (0_IK < bin .and. bin < size(timeZone%zone, kind = IK)) then
            abbr = trim(timeZone%Abbr(bin))
        elseif (zone == int(14 * 60, kind(zone))) then
            abbr = trim(timeZone%Abbr(size(timeZone%zone, kind = IK)))
        else
            abbr = SKG_""
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getHour12C
        use pm_kind, only: IKG => IK
        hour12 = getHour12(getHour())
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getHour12H
        use pm_kind, only: IKG => IK
        CHECK_ASSERTION(__LINE__, 0_IKG <= hour .and. hour <= 24_IKG, SK_"@getHour12(): The condition `0_IKG <= hour .and. hour <= 24_IKG` must hold. hour = "//getStr(hour)) ! fpp
        if (12_IKG < hour) then
            hour12 = hour - 12_IKG
        elseif (0_IKG < hour) then
            hour12 = hour
        else
            hour12 = 12_IKG
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getCountDaysInYear
        use pm_kind, only: IKG => IK
        if (isLeapYear(year)) then
            countDays = 366_IKG
        else
            countDays = 365_IKG
        end if
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getCountDaysInMonth
        use pm_kind, only: IKG => IK
        CHECK_ASSERTION(__LINE__, isValidDateTime(year, month), SK_"@getCountDays(): The input `year` and `month` must be a valid Gregorian date. year, month = "//getStr([year, month])) ! fpp
        if (isLeapYear(year)) then
            countDays = DAYS_OF_MONTH_LEAP(month)
        else
            countDays = DAYS_OF_MONTH(month)
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getCountWeeksInYear
        use pm_kind, only: IKG => IK
        integer(IKG) :: weekDayDec31Current, weekDayDec31Previous
        !weekDayDec31 = modulo(year + floor(0.25 * year, IKG) - floor(0.01 * year, IKG) + floor(0.0025 * year, IKG)), 7_IKG)
        weekDayDec31Current = getWeekDayISO(year, 12_IKG, 31_IKG)
        weekDayDec31Previous = getWeekDayISO(year - 1_IKG, 12_IKG, 31_IKG)
        if (weekDayDec31Current == 4_IKG .or. weekDayDec31Previous == 3_IKG) then ! current year ends Thursday or previous year ends Wednesday.
            countWeeks = 53_IKG
        else
            countWeeks = 52_IKG
        end if
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getCountWeeksInMonth
        use pm_kind, only: IKG => IK
        integer(IKG) :: countDays, weekdayStartISO
        CHECK_ASSERTION(__LINE__, isValidDateTime(year, month), SK_"@getCountWeeks(): The input `year` and `month` must be a valid Gregorian date. year, month = "//getStr([year, month])) ! fpp
        weekdayStartISO = getWeekDayISO(year, month, 1_IKG)
        countDays = getCountDays(year, month)
        if (countDays == 31_IKG) then
            if (weekdayStartISO < 2_IKG .or. weekdayStartISO > 4_IKG) then
                countWeeks = 4_IKG
            else
                countWeeks = 5_IKG
            end if
        elseif (countDays == 30_IKG) then
            if (weekdayStartISO < 3_IKG .or. weekdayStartISO > 4_IKG) then
                countWeeks = 4_IKG
            else
                countWeeks = 5_IKG
            end if
        elseif (countDays == 28_IKG) then
            countWeeks = 4_IKG
        elseif (countDays == 29_IKG) then
            if (weekdayStartISO /= 4_IKG) then
                countWeeks = 4_IKG
            else
                countWeeks = 5_IKG
            end if
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getWeekDayCurrent
        integer(IK) :: Values(8)
        call date_and_time(values = Values)
        CHECK_ASSERTION(__LINE__, all(Values /= -huge(0_IK)), SK_"@getWeekDay(): The processor does not possess a clock.")
        weekday = getWeekDay(year = Values(1), month = Values(2), day = Values(3))
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getWeekDayValues
        CHECK_ASSERTION(__LINE__, size(Values) > 2, SK_"@getWeekDay(): The condition `size(Values) > 2` must hold. size(Values) = "//getStr(size(Values)))
        CHECK_ASSERTION(__LINE__, size(Values) < 9, SK_"@getWeekDay(): The condition `size(Values) < 9` must hold. size(Values) = "//getStr(size(Values)))
        weekday = getWeekDay(year = Values(1), month = Values(2), day = Values(3))
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getWeekDayTriple
        use pm_kind, only: RKG => RK, IKG => IK
        !integer(IK) :: century
        !integer(IK) :: month_
        !integer(IK) :: year_
        CHECK_ASSERTION(__LINE__, isValidDateTime(year, month, day), SK_"@getWeekDay(): The input [year, month, day] must correspond to a valid Gregorian date. year, month, day = "//getStr([year, month, day]))
        weekday = modulo(floor(getJulianDay(year, month, day) + 1.5_RKG, IKG), 7_IKG)
        !if (month > 2_IK) then
        !    month_ = month + 1_IK
        !    year_ = year
        !else
        !    month_ = month + 13_IK
        !    year_ = year - 1_IK
        !end if
        !century = year_ / 100_IK
        !year_ = year_ - century * 100_IK
        !weekday = modulo(day + int(2.6 * real(month_), IK) + year_ + year_ / 4_IK + century / 4_IK - 2_IK * century, 7_IK)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getWeekYearCurrent
        use pm_kind, only: IKG => IK
        integer(IKG) :: WeekDate(3)
        WeekDate(1:3) = getWeekDate()
        weekYear = WeekDate(1)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getWeekYearValues
        use pm_kind, only: IKG => IK
        integer(IKG) :: WeekDate(3)
        WeekDate(1:3) = getWeekDate(Values)
        weekYear = WeekDate(1)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getWeekYearTriple
        use pm_kind, only: IKG => IK
        integer(IKG) :: WeekDate(3)
        WeekDate(1:3) = getWeekDate(year, month, day)
        weekYear = WeekDate(1)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getWeekDateCurrent
        integer(IK) :: Values(8)
        call date_and_time(values = Values)
        CHECK_ASSERTION(__LINE__, all(Values /= -huge(0_IK)), SK_"@getWeekDate(): The processor does not possess a clock.")
        WeekDate(1:3) = getWeekDate(year = Values(1), month = Values(2), day = Values(3))
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getWeekDateValues
        CHECK_ASSERTION(__LINE__, size(Values) > 2, SK_"@getWeekDate(): The condition `size(Values) > 2` must hold. size(Values) = "//getStr(size(Values)))
        CHECK_ASSERTION(__LINE__, size(Values) < 9, SK_"@getWeekDate(): The condition `size(Values) < 9` must hold. size(Values) = "//getStr(size(Values)))
        WeekDate(1:3) = getWeekDate(year = Values(1), month = Values(2), day = Values(3))
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getWeekDateTriple
        use pm_kind, only: IKG => IK
        CHECK_ASSERTION(__LINE__, isValidDateTime(year, month, day), SK_"@getWeekDate(): The input [year, month, day] must correspond to a valid Gregorian date. year, month, day = "//getStr([year, month, day]))
        WeekDate(1) = year
        WeekDate(3) = getWeekDayISO(year, month, day)
        WeekDate(2) = (10_IK + getOrdinalDay(year, month, day) - WeekDate(3)) / 7_IK
        if (WeekDate(2) == 0_IK) then
            WeekDate(1) = year - 1_IK
            WeekDate(2) = getCountWeeks(WeekDate(1))
        elseif (WeekDate(2) > getCountWeeks(year)) then
            WeekDate(1) = year + 1_IK
            WeekDate(2) = 1_IK
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getWeekDayISOCurrent
        integer(IK) :: Values(8)
        call date_and_time(values = Values)
        weekday = getWeekDayISO(year = Values(1), month = Values(2), day = Values(3))
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getWeekDayISOValues
        weekday = getWeekDayISO(year = Values(1), month = Values(2), day = Values(3))
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getWeekDayISOTriple
        weekday = getWeekDay(year, month, day)
        if (weekday == 0_IK) weekday = 7_IK
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getCountLeapYears
        integer(IK) :: year
        !check_assertion(__LINE__, until /= 0_IK, SK_"@getCountLeapYears(): The condition `until /= 0` must hold. until = "//getStr(until))
        if (until > 0_IK) then
            year = until
        else
            year = until - 4_IK
        end if
        countLeapYear = year / 4_IK - year / 100_IK + year / 400_IK
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getCountLeapYearsSince
        !check_assertion(__LINE__, until /= 0_IK, SK_"@getCountLeapYearsSince(): The condition `until /= 0` must hold. until = "//getStr(until))
        !check_assertion(__LINE__, since /= 0_IK, SK_"@getCountLeapYearsSince(): The condition `since /= 0` must hold. since = "//getStr(since))
        CHECK_ASSERTION(__LINE__, since <= until, SK_"@getCountLeapYears(): The condition `since <= until` must hold. since, until = "//getStr([since, until]))
        countLeapYear = getCountLeapYears(until) - getCountLeapYears(since)
        if (since > 0_IK .and. until > 0_IK) then
            if (isLeapYear(since)) countLeapYear = countLeapYear + 1_IK
        elseif (since <= 0_IK .and. until <= 0_IK) then
            if (isLeapYear(until)) countLeapYear = countLeapYear + 1_IK
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isLastDayInMonthC
        lastDayInMonth = isLastDayInMonth(getDateTime())
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isLastDayInMonthValues
        CHECK_ASSERTION(__LINE__, size(Values) > 2, SK_"@isLastDayInMonth(): The condition `size(Values) > 2` must hold. size(Values) = "//getStr(size(Values)))
        CHECK_ASSERTION(__LINE__, size(Values) < 9, SK_"@isLastDayInMonth(): The condition `size(Values) < 9` must hold. size(Values) = "//getStr(size(Values)))
        lastDayInMonth = isLastDayInMonth(year = Values(1), month = Values(2), day = Values(3))
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isLastDayInMonthTriple
        CHECK_ASSERTION(__LINE__, isValidDateTime(year, month, day), SK_"@isLastDayInMonth(): The input [year, month, day] must correspond to a valid Gregorian date. year, month, day = "//getStr([year, month, day]))
        if (month /= 2_IK) then
            lastDayInMonth = logical((month ==  1_IK .and. day == 31_IK) & ! LCOV_EXCL_LINE
                                .or. (month ==  3_IK .and. day == 31_IK) & ! LCOV_EXCL_LINE
                                .or. (month ==  4_IK .and. day == 30_IK) & ! LCOV_EXCL_LINE
                                .or. (month ==  5_IK .and. day == 31_IK) & ! LCOV_EXCL_LINE
                                .or. (month ==  6_IK .and. day == 30_IK) & ! LCOV_EXCL_LINE
                                .or. (month ==  7_IK .and. day == 31_IK) & ! LCOV_EXCL_LINE
                                .or. (month ==  8_IK .and. day == 31_IK) & ! LCOV_EXCL_LINE
                                .or. (month ==  9_IK .and. day == 30_IK) & ! LCOV_EXCL_LINE
                                .or. (month == 10_IK .and. day == 31_IK) & ! LCOV_EXCL_LINE
                                .or. (month == 11_IK .and. day == 30_IK) & ! LCOV_EXCL_LINE
                                .or. (month == 12_IK .and. day == 31_IK), LK)
        else
            if (isLeapYear(year)) then
                lastDayInMonth = logical(day == 29_IK, LK)
            else
                lastDayInMonth = logical(day == 28_IK, LK)
            end if
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getDateAfterC
        DateAfter(1:3) = getDateAfter(getDateTime())
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getDateAfterValues
        CHECK_ASSERTION(__LINE__, size(Values) > 2, SK_"@getDateAfter(): The condition `size(Values) > 2` must hold. size(Values) = "//getStr(size(Values)))
        CHECK_ASSERTION(__LINE__, size(Values) < 9, SK_"@getDateAfter(): The condition `size(Values) < 9` must hold. size(Values) = "//getStr(size(Values)))
        DateAfter(1:3) = getDateAfter(year = Values(1), month = Values(2), day = Values(3))
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getDateAfterTriple
        use pm_kind, only: IKG => IK
        CHECK_ASSERTION(__LINE__, isValidDateTime(year, month, day), SK_"@getDateAfter(): The input [year, month, day] must correspond to a valid Gregorian date. year, month, day = "//getStr([year, month, day]))
        if (isLastDayInMonth(year, month, day)) then
            if (month < 12_IK) then
                DateAfter(1) = year
                DateAfter(2) = month + 1_IK
                DateAfter(3) = 1_IK
            else
                DateAfter(1) = year + 1_IKG
                !if (year /= -1_IKG) then
                !    DateAfter(1) = year + 1_IKG
                !else
                !    DateAfter(1) = 1_IKG
                !end if
                DateAfter(2) = 1_IK
                DateAfter(3) = 1_IK
            end if
        else
            DateAfter(1) = year
            DateAfter(2) = month
            DateAfter(3) = day + 1_IK
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getDateBeforeC
        DateBefore(1:3) = getDateBefore(getDateTime())
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getDateBeforeValues
        CHECK_ASSERTION(__LINE__, size(Values) > 2, SK_"@getDateBefore(): The condition `size(Values) > 2` must hold. size(Values) = "//getStr(size(Values)))
        CHECK_ASSERTION(__LINE__, size(Values) < 9, SK_"@getDateBefore(): The condition `size(Values) < 9` must hold. size(Values) = "//getStr(size(Values)))
        DateBefore(1:3) = getDateBefore(year = Values(1), month = Values(2), day = Values(3))
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getDateBeforeTriple
        use pm_kind, only: IKG => IK
        CHECK_ASSERTION(__LINE__, isValidDateTime(year, month, day), SK_"@getDateBefore(): The input [year, month, day] must correspond to a valid Gregorian date. year, month, day = "//getStr([year, month, day]))
        if (day > 1_IKG) then
            DateBefore(1) = year
            DateBefore(2) = month
            DateBefore(3) = day - 1_IKG
        else
            if (month > 1_IKG) then
                DateBefore(1) = year
                DateBefore(2) = month - 1_IKG
                if (isLeapYear(year)) then
                    DateBefore(3) = DAYS_OF_MONTH_LEAP(DateBefore(2))
                else
                    DateBefore(3) = DAYS_OF_MONTH(DateBefore(2))
                end if
            else ! month == 1_IKG
                !if (year /= 1_IKG) then
                !    DateBefore(1) = year - 1_IKG
                !else
                !    DateBefore(1) = -1_IKG ! leap year
                !end if
                DateBefore(1) = year - 1_IKG
                DateBefore(2) = 12_IKG
                DateBefore(3) = 31_IKG ! December always has 31 days.
            end if
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getOrdinalDayCurrent
        integer(IK) :: Values(8)
        call date_and_time(values = Values)
        ordinalDay = getOrdinalDay(year = Values(1), month = Values(2), day = Values(3))
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getOrdinalDayValues
        CHECK_ASSERTION(__LINE__, size(Values) > 2, SK_"@getOrdinalDay(): The condition `size(Values) > 2` must hold. size(Values) = "//getStr(size(Values)))
        CHECK_ASSERTION(__LINE__, size(Values) < 9, SK_"@getOrdinalDay(): The condition `size(Values) < 9` must hold. size(Values) = "//getStr(size(Values)))
        ordinalDay = getOrdinalDay(year = Values(1), month = Values(2), day = Values(3))
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getOrdinalDayTriple
        use pm_kind, only: IKG => IK
        integer(IKG), parameter :: CUMULATIVE_MONTH_DAYS(12) = [0_IKG, 31_IKG, 59_IKG, 90_IKG, 120_IKG, 151_IKG, 181_IKG, 212_IKG, 243_IKG, 273_IKG, 304_IKG, 334_IKG]
        CHECK_ASSERTION(__LINE__, isValidDateTime(year, month, day), SK_"@getOrdinalDay(): The input [year, month, day] must correspond to a valid Gregorian date. year, month, day = "//getStr([year, month, day]))
        ordinalDay = CUMULATIVE_MONTH_DAYS(month) + day
        if (isLeapYear(year) .and. month > 2_IK) ordinalDay = ordinalDay + 1_IK
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getWeekNumberCurrent
        integer(IK) :: Values(8)
        call date_and_time(values = Values)
        weekNumber = getWeekNumber(year = Values(1), month = Values(2), day = Values(3))
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getWeekNumberValues
        CHECK_ASSERTION(__LINE__, size(Values) > 2, SK_"@getWeekNumber(): The condition `size(Values) > 2` must hold. size(Values) = "//getStr(size(Values)))
        CHECK_ASSERTION(__LINE__, size(Values) < 9, SK_"@getWeekNumber(): The condition `size(Values) < 9` must hold. size(Values) = "//getStr(size(Values)))
        weekNumber = getWeekNumber(year = Values(1), month = Values(2), day = Values(3))
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getWeekNumberTriple
        use pm_kind, only: IKG => IK
        integer(IKG), parameter :: CUMULATIVE_MONTH_DAYS(12) = [0_IKG, 31_IKG, 59_IKG, 90_IKG, 120_IKG, 151_IKG, 181_IKG, 212_IKG, 243_IKG, 273_IKG, 304_IKG, 334_IKG]
        CHECK_ASSERTION(__LINE__, isValidDateTime(year, month, day), SK_"@getWeekNumber(): The input [year, month, day] must correspond to a valid Gregorian date. year, month, day = "//getStr([year, month, day]))
        weekNumber = (10_IK + getOrdinalDay(year, month, day) - getWeekDayISO(year, month, day)) / 7_IK
        if (weekNumber == 0_IK) then
            weekNumber = getCountWeeks(year - 1_IK)
        elseif (weekNumber > getCountWeeks(year)) then
            weekNumber = 1_IK
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getDateTimeDiffValues
        CHECK_ASSERTION(__LINE__, size(Values1) == 8, SK_"@getDateTimeDiff(): The condition `size(Values1) == 8` must hold. size(Values1) = "//getStr(size(Values1)))
        CHECK_ASSERTION(__LINE__, size(Values2) == 8, SK_"@getDateTimeDiff(): The condition `size(Values2) == 8` must hold. size(Values2) = "//getStr(size(Values2)))
        if (Values1(4) == Values2(4)) then ! identical timezone.
            dateTimeDiff= getJulianDay(Values1(1), Values1(2), Values1(3), 0_IK, Values1(5), Values1(6), Values1(7), Values1(8)) & ! LCOV_EXCL_LINE
                        - getJulianDay(Values2(1), Values2(2), Values2(3), 0_IK, Values2(5), Values2(6), Values2(7), Values2(8))
        else ! different timezones.
            dateTimeDiff= getJulianDay(Values1(1), Values1(2), Values1(3), Values1(4), Values1(5), Values1(6), Values1(7), Values1(8)) & ! LCOV_EXCL_LINE
                        - getJulianDay(Values2(1), Values2(2), Values2(3), Values2(4), Values2(5), Values2(6), Values2(7), Values2(8))
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isMorningCurrent
        use pm_kind, only: IKG => IK
        integer(IKG) :: Values(8)
        call date_and_time(values = Values)
        morning = Values(5) < 12_IKG
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isMorningZ
        use pm_kind, only: IKG => IK
        morning = getHour(zone) < 12_IKG
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isMorningJD
        use pm_kind, only: RKG => RK
        morning = logical(julianDay - real(floor(julianDay, IK), RKG) >= 0.5_RKG, LK)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isMorningJDZ
        use pm_kind, only: RKG => RK
        real(RKG), parameter :: INV_MINUTES_PER_DAY = 1._RKG / MINUTES_PER_DAY
        CHECK_ASSERTION(__LINE__, isValidZone(zone), SK_"@isMorning(): The condition `isValidZone(zone)` must hold. zone = "//getStr(zone)) ! fpp
        morning = isMorning(julianDay + zone * INV_MINUTES_PER_DAY)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getJulianDay_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getJulianDayC_ENABLED 1

    module procedure getJulianDayC
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getJulianDayC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getJulianDayV_ENABLED 1

    module procedure getJulianDayV
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getJulianDayV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getJulianDayY_ENABLED 1

    module procedure getJulianDayY
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getJulianDayY_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getJulianDayYM_ENABLED 1

    module procedure getJulianDayYM
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getJulianDayYM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getJulianDayYMD_ENABLED 1

    module procedure getJulianDayYMD
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getJulianDayYMD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getJulianDayYMDZ_ENABLED 1

    module procedure getJulianDayYMDZ
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getJulianDayYMDZ_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getJulianDayYMDZH_ENABLED 1

    module procedure getJulianDayYMDZH
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getJulianDayYMDZH_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getJulianDayYMDZHM_ENABLED 1

    module procedure getJulianDayYMDZHM
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getJulianDayYMDZHM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getJulianDayYMDZHMS_ENABLED 1

    module procedure getJulianDayYMDZHMS
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getJulianDayYMDZHMS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getJulianDayYMDZHMSM_ENABLED 1

    module procedure getJulianDayYMDZHMSM
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getJulianDayYMDZHMSM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getJulianDay_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeUTC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeUTCC_ENABLED 1

    module procedure getDateTimeUTCC
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeUTCC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeUTCV_ENABLED 1

    module procedure getDateTimeUTCV
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeUTCV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeUTCYMDZ_ENABLED 1

    module procedure getDateTimeUTCYMDZ
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeUTCYMDZ_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeUTCYMDZH_ENABLED 1

    module procedure getDateTimeUTCYMDZH
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeUTCYMDZH_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeUTCYMDZHM_ENABLED 1

    module procedure getDateTimeUTCYMDZHM
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeUTCYMDZHM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeUTCYMDZHMS_ENABLED 1

    module procedure getDateTimeUTCYMDZHMS
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeUTCYMDZHMS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeUTCYMDZHMSM_ENABLED 1

    module procedure getDateTimeUTCYMDZHMSM
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeUTCYMDZHMSM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDateTimeUTC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeNewZone_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeNewZoneC_ENABLED 1

    module procedure getDateTimeNewZoneC
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeNewZoneC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeNewZoneV_ENABLED 1

    module procedure getDateTimeNewZoneV
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeNewZoneV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeNewZoneYMDZ_ENABLED 1

    module procedure getDateTimeNewZoneYMDZ
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeNewZoneYMDZ_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeNewZoneYMDZH_ENABLED 1

    module procedure getDateTimeNewZoneYMDZH
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeNewZoneYMDZH_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeNewZoneYMDZHM_ENABLED 1

    module procedure getDateTimeNewZoneYMDZHM
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeNewZoneYMDZHM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeNewZoneYMDZHMS_ENABLED 1

    module procedure getDateTimeNewZoneYMDZHMS
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeNewZoneYMDZHMS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeNewZoneYMDZHMSM_ENABLED 1

    module procedure getDateTimeNewZoneYMDZHMSM
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeNewZoneYMDZHMSM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDateTimeNewZone_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeShifted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeShiftedC_ENABLED 1

    module procedure getDateTimeShiftedC
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeShiftedC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeShiftedV_ENABLED 1

    module procedure getDateTimeShiftedV
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeShiftedV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeShiftedY_ENABLED 1

    module procedure getDateTimeShiftedY
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeShiftedY_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeShiftedYM_ENABLED 1

    module procedure getDateTimeShiftedYM
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeShiftedYM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeShiftedYMD_ENABLED 1

    module procedure getDateTimeShiftedYMD
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeShiftedYMD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeShiftedYMDZ_ENABLED 1

    module procedure getDateTimeShiftedYMDZ
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeShiftedYMDZ_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeShiftedYMDZH_ENABLED 1

    module procedure getDateTimeShiftedYMDZH
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeShiftedYMDZH_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeShiftedYMDZHM_ENABLED 1

    module procedure getDateTimeShiftedYMDZHM
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeShiftedYMDZHM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeShiftedYMDZHMS_ENABLED 1

    module procedure getDateTimeShiftedYMDZHMS
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeShiftedYMDZHMS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeShiftedYMDZHMSM_ENABLED 1

    module procedure getDateTimeShiftedYMDZHMSM
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeShiftedYMDZHMSM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDateTimeShifted_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isValidDateTime_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isValidDateTimeV_ENABLED 1

    module procedure isValidDateTimeV
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef isValidDateTimeV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isValidDateTimeY_ENABLED 1

    module procedure isValidDateTimeY
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef isValidDateTimeY_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isValidDateTimeYM_ENABLED 1

    module procedure isValidDateTimeYM
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef isValidDateTimeYM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isValidDateTimeYMD_ENABLED 1

    module procedure isValidDateTimeYMD
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef isValidDateTimeYMD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isValidDateTimeYMDZ_ENABLED 1

    module procedure isValidDateTimeYMDZ
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef isValidDateTimeYMDZ_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isValidDateTimeYMDZH_ENABLED 1

    module procedure isValidDateTimeYMDZH
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef isValidDateTimeYMDZH_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isValidDateTimeYMDZHM_ENABLED 1

    module procedure isValidDateTimeYMDZHM
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef isValidDateTimeYMDZHM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isValidDateTimeYMDZHMS_ENABLED 1

    module procedure isValidDateTimeYMDZHMS
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef isValidDateTimeYMDZHMS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isValidDateTimeYMDZHMSM_ENABLED 1

    module procedure isValidDateTimeYMDZHMSM
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef isValidDateTimeYMDZHMSM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isValidDateTime_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTime_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeValues_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeValuesJ_ENABLED 1

    module procedure getDateTimeValuesJ
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeValuesJ_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeValuesJZ_ENABLED 1

    module procedure getDateTimeValuesJZ
        use pm_kind, only: IKG => IK, RKG => RK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeValuesJZ_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeValuesC_ENABLED 1

    module procedure getDateTimeValuesC
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeValuesC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeValuesY_ENABLED 1

    module procedure getDateTimeValuesY
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeValuesY_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeValuesYM_ENABLED 1

    module procedure getDateTimeValuesYM
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeValuesYM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeValuesYMD_ENABLED 1

    module procedure getDateTimeValuesYMD
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeValuesYMD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeValuesYMDZ_ENABLED 1

    module procedure getDateTimeValuesYMDZ
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeValuesYMDZ_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeValuesYMDZH_ENABLED 1

    module procedure getDateTimeValuesYMDZH
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeValuesYMDZH_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeValuesYMDZHM_ENABLED 1

    module procedure getDateTimeValuesYMDZHM
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeValuesYMDZHM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeValuesYMDZHMS_ENABLED 1

    module procedure getDateTimeValuesYMDZHMS
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeValuesYMDZHMS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeValuesYMDZHMSM_ENABLED 1

    module procedure getDateTimeValuesYMDZHMSM
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeValuesYMDZHMSM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDateTimeValues_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeString_ENABLED 1

!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getDateTimeStringJ_ENABLED 1
!
!    module procedure getDateTimeStringJ
!        use pm_kind, only: IKG => IK, RKG => RK
!#include "pm_dateTime@routines.inc.F90"
!    end procedure
!
!#undef getDateTimeStringJ_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getDateTimeStringJZ_ENABLED 1
!
!    module procedure getDateTimeStringJZ
!        use pm_kind, only: IKG => IK, RKG => RK
!#include "pm_dateTime@routines.inc.F90"
!    end procedure
!
!#undef getDateTimeStringJZ_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeStringC_ENABLED 1

    module procedure getDateTimeStringC
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeStringC_ENABLED

!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDateTimeStringV_ENABLED 1

    module procedure getDateTimeStringV
        use pm_kind, only: IKG => IK, SKG => SK
#include "pm_dateTime@routines.inc.F90"
    end procedure

#undef getDateTimeStringV_ENABLED

!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getDateTimeStringY_ENABLED 1
!
!    module procedure getDateTimeStringY
!        use pm_kind, only: IKG => IK, SKG => SK
!#include "pm_dateTime@routines.inc.F90"
!    end procedure
!
!#undef getDateTimeStringY_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getDateTimeStringYM_ENABLED 1
!
!    module procedure getDateTimeStringYM
!        use pm_kind, only: IKG => IK, SKG => SK
!#include "pm_dateTime@routines.inc.F90"
!    end procedure
!
!#undef getDateTimeStringYM_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getDateTimeStringYMD_ENABLED 1
!
!    module procedure getDateTimeStringYMD
!        use pm_kind, only: IKG => IK, SKG => SK
!#include "pm_dateTime@routines.inc.F90"
!    end procedure
!
!#undef getDateTimeStringYMD_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getDateTimeStringYMDZ_ENABLED 1
!
!    module procedure getDateTimeStringYMDZ
!        use pm_kind, only: IKG => IK, SKG => SK
!#include "pm_dateTime@routines.inc.F90"
!    end procedure
!
!#undef getDateTimeStringYMDZ_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getDateTimeStringYMDZH_ENABLED 1
!
!    module procedure getDateTimeStringYMDZH
!        use pm_kind, only: IKG => IK, SKG => SK
!#include "pm_dateTime@routines.inc.F90"
!    end procedure
!
!#undef getDateTimeStringYMDZH_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getDateTimeStringYMDZHM_ENABLED 1
!
!    module procedure getDateTimeStringYMDZHM
!        use pm_kind, only: IKG => IK, SKG => SK
!#include "pm_dateTime@routines.inc.F90"
!    end procedure
!
!#undef getDateTimeStringYMDZHM_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getDateTimeStringYMDZHMS_ENABLED 1
!
!    module procedure getDateTimeStringYMDZHMS
!        use pm_kind, only: IKG => IK, SKG => SK
!#include "pm_dateTime@routines.inc.F90"
!    end procedure
!
!#undef getDateTimeStringYMDZHMS_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getDateTimeStringYMDZHMSM_ENABLED 1
!
!    module procedure getDateTimeStringYMDZHMSM
!        use pm_kind, only: IKG => IK, SKG => SK
!#include "pm_dateTime@routines.inc.F90"
!    end procedure
!
!#undef getDateTimeStringYMDZHMSM_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDateTimeString_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDateTime_ENABLED

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines

#undef  CHECK_ASSERTION