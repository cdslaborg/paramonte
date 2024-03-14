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
!>  This include file contains procedure implementations of [pm_dateTime](@ref pm_dateTime).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, March 23, 2012, 5:33 PM, National Institute for Fusion Studies, The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%
#if     isValidZone_ENABLED
        !%%%%%%%%%%%%%%%%%%

        isValid = logical(int(ZONE_MIN, IKC) <= zone .and. zone <= int(ZONE_MAX, IKC), LK)

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   getDateTimeDiff_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, size(Values1, kind = IK) == 8_IK, SK_"@getDateTimeDiff(): The condition `size(Values1) == 8` must hold. size(Values1) = "//getStr(size(Values1, kind = IK)))
        CHECK_ASSERTION(__LINE__, size(Values2, kind = IK) == 8_IK, SK_"@getDateTimeDiff(): The condition `size(Values2) == 8` must hold. size(Values2) = "//getStr(size(Values2, kind = IK)))
        if (Values1(4) == Values2(4)) then ! identical timezone.
            dateTimeDiff= getJulianDay(Values1(1), Values1(2), Values1(3), 0_IK, Values1(5), Values1(6), Values1(7), Values1(8)) & ! LCOV_EXCL_LINE
                        - getJulianDay(Values2(1), Values2(2), Values2(3), 0_IK, Values2(5), Values2(6), Values2(7), Values2(8))
        else ! different timezones.
            dateTimeDiff= getJulianDay(Values1(1), Values1(2), Values1(3), Values1(4), Values1(5), Values1(6), Values1(7), Values1(8)) & ! LCOV_EXCL_LINE
                        - getJulianDay(Values2(1), Values2(2), Values2(3), Values2(4), Values2(5), Values2(6), Values2(7), Values2(8))
        end if

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   isValidDateTime_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

#if     isValidDateTimeV_ENABLED
        integer(IK) :: lenValues
        lenValues = size(values, 1, IK)
        if (lenValues == 8_IK) then
            isValid = isValidDateTime(values(1), values(2), values(3), values(4), values(5), values(6), values(7), values(8))
        elseif (lenValues == 7_IK) then
            isValid = isValidDateTime(values(1), values(2), values(3), values(4), values(5), values(6), values(7))
        elseif (lenValues == 6_IK) then
            isValid = isValidDateTime(values(1), values(2), values(3), values(4), values(5), values(6))
        elseif (lenValues == 5_IK) then
            isValid = isValidDateTime(values(1), values(2), values(3), values(4), values(5))
        elseif (lenValues == 4_IK) then
            isValid = isValidDateTime(values(1), values(2), values(3), values(4))
        elseif (lenValues == 3_IK) then
            isValid = isValidDateTime(values(1), values(2), values(3))
        elseif (lenValues == 2_IK) then
            isValid = isValidDateTime(values(1), values(2))
        elseif (lenValues == 1_IK) then
            isValid = isValidDateTime(values(1))
        else
            isValid = .false._LK
            !error stop MODULE_NAME//SK_"@isValidDateTime(): The length of the input argument `values` must be less than 9 and non-zero. size(values) = "//getStr(lenValues)
        end if
#elif   isValidDateTimeY_ENABLED
        isValid = .true._LK ! year /= 0_IKC ! year zero is explicitly allowed in ISO 8601.
#elif   isValidDateTimeYM_ENABLED
        isValid = isValidDateTime(year) .and. 0_IKC < month .and. month < 13_IKC
#elif   isValidDateTimeYMD_ENABLED
        isValid = isValidDateTime(year, month) .and. day > 0_IKC
        if (isValid) then
            if (isLeapYear(year)) then
                isValid = day <= DAYS_OF_MONTH_LEAP(month)
            else
                isValid = day <= DAYS_OF_MONTH(month)
            end if
        end if
#elif   isValidDateTimeYMDZ_ENABLED
        isValid = isValidDateTime(year, month, day) .and. isValidZone(zone)
#elif   isValidDateTimeYMDZH_ENABLED
        isValid = isValidDateTime(year, month, day, zone) .and. 0_IKC <= hour .and. hour < 24_IKC
#elif   isValidDateTimeYMDZHM_ENABLED
        isValid = isValidDateTime(year, month, day, zone, hour) .and. 0_IKC <= minute .and. minute < 60_IKC
#elif   isValidDateTimeYMDZHMS_ENABLED
        isValid = isValidDateTime(year, month, day, zone, hour, minute) .and. 0_IKC <= second .and. second < 60_IKC
#elif   isValidDateTimeYMDZHMSM_ENABLED
        isValid = isValidDateTime(year, month, day, zone, hour, minute, second) .and. 0_IKC <= millisecond .and. millisecond < 1000_IKC
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getDateTimeValues_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

#if     getDateTimeValuesJ_ENABLED
        real(RKC)   , parameter :: K1 = 0.25_RKC ! The condition `0.002929687499688476 < K1 <= 0.2521972656249999` must hold.
        real(RKC)   , parameter :: K2 = 0.25_RKC ! The condition `0. < K2 <= 0.25` must hold.
        real(RKC)   , parameter :: MEAN_YEAR = 365.25_RKC
        real(RKC)   , parameter :: MEAN_YEAR_INVERSE = 1._RKC / MEAN_YEAR
        real(RKC)   , parameter :: WHOLE_CENTURY_FACTOR = 1._RKC / 36524.25_RKC
        integer(IKC), parameter :: DAYS_IN_PAST_MONTHS(3:14) = [0_IKC, 31_IKC, 61_IKC, 92_IKC, 122_IKC, 153_IKC, 184_IKC, 214_IKC, 245_IKC, 275_IKC, 306_IKC, 337_IKC] ! last two correspond to jan, feb.
        integer(IKC)    :: julianDayOffset_IK ! Z
        integer(IKC)    :: fullCenturyCount
        integer(IKC)    :: dayOfYearTerm
        real(RKC)       :: hours_RK
        real(RKC)       :: minutes_RK
        real(RKC)       :: seconds_RK
        real(RKC)       :: julianDayOffset_RK
        real(RKC)       :: julianDayOffsetResidual ! R
        real(RKC)       :: daysInWholeCenturyMinusConst ! B
        julianDayOffset_RK = julianDay - 1721118.5_RKC
        julianDayOffset_IK = floor(julianDayOffset_RK, IKC) ! Z
        julianDayOffsetResidual = julianDayOffset_RK - real(julianDayOffset_IK, RKC) ! R : always positive
        fullCenturyCount = floor((real(julianDayOffset_IK, RKC) - K1) * WHOLE_CENTURY_FACTOR, IKC) ! A
        dayOfYearTerm = julianDayOffset_IK + fullCenturyCount - floor(fullCenturyCount * 0.25_RKC, IKC) ! Z + A - floor(A/4)
        daysInWholeCenturyMinusConst = dayOfYearTerm - K2 ! B
        values(1) = floor(daysInWholeCenturyMinusConst * MEAN_YEAR_INVERSE, IKC) ! Y : Calendar Year Starting March.
        values(3) = dayOfYearTerm - floor(values(1) * MEAN_YEAR, IKC) ! C : Day Of Year.
        values(2) = (5_IKC * values(3) + 456_IKC) / 153_IKC ! M : Month Of Year in the range 3:14.
        values(3) = values(3) - DAYS_IN_PAST_MONTHS(values(2)) !+ int(julianDayOffsetResidual, IKC)
        !if (julianDayOffsetResidual > 1._RKC) error stop getStr(julianDayOffsetResidual)
        if (values(2) > 12_IKC) then
            values(1) = values(1) + 1_IKC
            values(2) = values(2) - 12_IKC
        end if
        values(4) = 0_IKC
        hours_RK = julianDayOffsetResidual * 24._RKC
        values(5) = int(hours_RK, IKC)
        !if (values(5) < 0_IKC) then
        !    values(1:3) = getDateBefore(values(1), values(2), values(3))
        !    values(5) = 24_IKC - values(5)
        !end if
        minutes_RK = (hours_RK - values(5)) * 60._RKC
        values(6) = int(minutes_RK, IKC)
        seconds_RK = (minutes_RK - values(6)) * 60._RKC
        values(7) = int(seconds_RK, IKC)
        values(8) = nint((seconds_RK - values(7)) * 1000._RKC, IKC)
        if (values(8) == 1000_IKC) then
            values(8) = 0_IKC
            values(7) = values(7) + 1_IKC
            if (values(7) == 60_IKC) then
                values(7) = 0_IKC
                values(6) = values(6) + 1_IKC
                if (values(6) == 60_IKC) then
                    values(6) = 0_IKC
                    values(5) = values(5) + 1_IKC
                end if
                if (values(5) == 24_IKC) then
                    values(5) = 0_IKC
                    values(1:3) = getDateAfter(values(1), values(2), values(3))
                end if
            end if
        end if
#elif   getDateTimeValuesJZ_ENABLED
        values(1:8) = getDateTimeNewZone(zone, getDateTime(julianDay))
#elif   getDateTimeValuesC_ENABLED
        call date_and_time(values = values)
        !integer(IK)     :: lenValues
        !integer(IKC)    :: Values_(8)
        !lenValues = size(values, 1, IK)
        !CHECK_ASSERTION(__LINE__, 0_IK < lenValues .and. lenValues < 9_IK, SK_"@getDateTime(): The input argument `values` must have a non-zero size that is less than 9. size(values) = "//getStr(lenValues)) ! fpp
        !if (lenValues > 7_IK) then
        !    call date_and_time(values = values)
        !else
        !    call date_and_time(values = Values_)
        !    values(1:lenValues) = Values_(1:lenValues)
        !end if
        !CHECK_ASSERTION(__LINE__, all(values(1:lenValues) /= -huge(0_IKC)), SK_"@getDateTime(): The processor does not have a clock.") ! fpp
#elif   getDateTimeValuesY_ENABLED
        values(1) = year
        values(2) = 1_IKC
        values(3) = 1_IKC
        values(4:8) = 0_IKC
#elif   getDateTimeValuesYM_ENABLED
        values(1) = year
        values(2) = month
        values(3) = 1_IKC
        values(4:8) = 0_IKC
#elif   getDateTimeValuesYMD_ENABLED
        values(1) = year
        values(2) = month
        values(3) = day
        values(4:8) = 0_IKC
#elif   getDateTimeValuesYMDZ_ENABLED
        values(1) = year
        values(2) = month
        values(3) = day
        values(4) = zone
        values(5:8) = 0_IKC
#elif   getDateTimeValuesYMDZH_ENABLED
        values(1) = year
        values(2) = month
        values(3) = day
        values(4) = zone
        values(5) = hour
        values(6:8) = 0_IKC
#elif   getDateTimeValuesYMDZHM_ENABLED
        values(1) = year
        values(2) = month
        values(3) = day
        values(4) = zone
        values(5) = hour
        values(6) = minute
        values(7:8) = 0_IKC
#elif   getDateTimeValuesYMDZHMS_ENABLED
        values(1) = year
        values(2) = month
        values(3) = day
        values(4) = zone
        values(5) = hour
        values(6) = minute
        values(7) = second
        values(8) = 0_IKC
#elif   getDateTimeValuesYMDZHMSM_ENABLED
        values(1) = year
        values(2) = month
        values(3) = day
        values(4) = zone
        values(5) = hour
        values(6) = minute
        values(7) = second
        values(8) = millisecond
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getDateTimeString_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

#define RESIZE_STRING(LENSEG) \
eposnew = epos + LENSEG; \
if (eposnew > lenString) then; \
if (allocated(tempstr)) deallocate(tempstr); \
allocate(character(eposnew,SKC) :: tempstr); \
tempstr(1:lenString) = string; \
call move_alloc(tempstr, string); \
lenString = eposnew; \
end if;

#if     getDateTimeStringC_ENABLED
        integer(IKC)    :: values(8)
        call date_and_time(values = values)
        string = getDateTime(format, values)
#elif   getDateTimeStringV_ENABLED
        !>  \warning
        !>  The output of getStr() in this procedure is of kind \SK which is incompatible with any value of SKC /= SK.
        !>  For now, this is not an issue since both kinds point to the default character kind.
        !>  This will however become an issue once Fortran standard and compilers support non-default date and time characters.
        use pm_val2str, only: getStr
        character(:,SKC), allocatable   :: tempstr, abbr
       !character(1,SKC)                :: sep
       !character(28,SKC)               :: workspace
        character(9,SKC)                :: workspace9
        integer(IK)                     :: lenString, lenFormat, i, epos, eposnew, lenSeg
        integer(IKC)                    :: century!, WeekDate(3)
        allocate(character(127,SKC)     :: string)
        lenFormat = len(format, IK)
        lenString = 127_IK
        eposnew = 0_IK ! the last touched (end) position in the string
        i = 0_IK
        CHECK_ASSERTION(__LINE__, 0 < size(values) .and. size(values) < 9, SK_"@getDateTime(): The condition `0 < size(values) .and. size(values) < 9` must hold. size(values) = "//getStr(size(values))) ! fpp
        CHECK_ASSERTION(__LINE__, merge(size(values) > 1, .true., \
        index(format, "%b") > 0 .or. \
        index(format, "%B") > 0 .or. \
        index(format, "%h") > 0 .or. \
        index(format, "%m") > 0 .or. \
        index(format, "%y") > 0 \
        ), SK_"@getDateTime(): The condition `size(values) > 1` must hold. size(values) = "//getStr(size(values))) ! fpp
        CHECK_ASSERTION(__LINE__, merge(size(values) > 2, .true., \
        index(format, "%a") > 0 .or. \
        index(format, "%A") > 0 .or. \
        index(format, "%d") > 0 .or. \
        index(format, "%D") > 0 .or. \
        index(format, "%e") > 0 .or. \
        index(format, "%g") > 0 .or. \
        index(format, "%G") > 0 .or. \
        index(format, "%j") > 0 .or. \
        index(format, "%u") > 0 .or. \
        index(format, "%U") > 0 .or. \
        index(format, "%V") > 0 .or. \
        index(format, "%w") > 0 .or. \
        index(format, "%W") > 0 .or. \
        index(format, "%x") > 0 \
        ), SK_"@getDateTime(): The condition `size(values) > 2` must hold. size(values) = "//getStr(size(values))) ! fpp
        CHECK_ASSERTION(__LINE__, merge(size(values) > 3, .true., \
        index(format, "%z") > 0 .or. \
        index(format, "%Z") > 0 \
        ), SK_"@getDateTime(): The condition `size(values) > 4` must hold. size(values) = "//getStr(size(values))) ! fpp
        CHECK_ASSERTION(__LINE__, merge(size(values) > 4, .true., \
        index(format, "%H") > 0 .or. \
        index(format, "%I") > 0 .or. \
        index(format, "%p") > 0 \
        ), SK_"@getDateTime(): The condition `size(values) > 4` must hold. size(values) = "//getStr(size(values))) ! fpp
        CHECK_ASSERTION(__LINE__, merge(size(values) > 5, .true., \
        index(format, "%M") > 0 .or. \
        index(format, "%R") > 0 \
        ), SK_"@getDateTime(): The condition `size(values) > 5` must hold. size(values) = "//getStr(size(values))) ! fpp
        CHECK_ASSERTION(__LINE__, merge(size(values) > 6, .true., \
        index(format, "%c") > 0 .or. \
        index(format, "%r") > 0 .or. \
        index(format, "%S") > 0 .or. \
        index(format, "%T") > 0 .or. \
        index(format, "%X") > 0 \
        ), SK_"@getDateTime(): The condition `size(values) > 6` must hold. size(values) = "//getStr(size(values))) ! fpp
        CHECK_ASSERTION(__LINE__, merge(size(values) > 7, .true., index(format, "%f") > 0), SK_"@getDateTime(): The condition `size(values) > 7` must hold. size(values) = "//getStr(size(values))) ! fpp
        do
            i = i + 1_IK
            if (i > lenFormat) exit
            epos = eposnew
            if (format(i:i) == SKC_"%") then
                if (i == lenFormat) exit
                i = i + 1_IK
                if (format(i:i) == SKC_"a") then ! Abbreviated weekday name *
                    lenSeg = 3_IK
                    RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = WEEKDAY_NAME_ISO(getWeekDayISO(values(1), values(2), values(3)))(1:lenSeg)
                elseif (format(i:i) == SKC_"A") then ! Full weekday name *
                    workspace9 = WEEKDAY_NAME_ISO(getWeekDayISO(values(1), values(2), values(3)))
                    lenSeg = len_trim(workspace9, IKC)
                    RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = workspace9(1:lenSeg)
                elseif (format(i:i) == SKC_"b" .or. format(i:i) == SKC_"h") then ! Abbreviated month name * .or. Abbreviated month name * (same as %b)
                    lenSeg = 3_IK
                    RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = MONTH_NAME(values(2))(1:lenSeg)
                elseif (format(i:i) == SKC_"B") then ! Full month name *
                    lenSeg = len_trim(MONTH_NAME(values(2)), IKC)
                    RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = MONTH_NAME(values(2))(1:lenSeg)
                elseif (format(i:i) == SKC_"c") then ! Date and time representation *
                    if (values(1) > 0_IKC) then
                        RESIZE_STRING(24_IK) ! fpp sets eposnew, and resizes string.
                    else
                        RESIZE_STRING(25_IK) ! fpp sets eposnew, and resizes string.
                    end if
                    string(epos + 1 : eposnew) =    WEEKDAY_NAME_ISO(getWeekDayISO(values(1), values(2), values(3)))(1:3)//SKC_" "// & ! LCOV_EXCL_LINE
                                                    MONTH_NAME(values(2))(1:3)//SKC_" "// & ! LCOV_EXCL_LINE
                                                    getStr(values(3), length = 2_IK, format = "(1I0.2)")//SKC_" "// & ! LCOV_EXCL_LINE
                                                    getStr(values(5), length = 2_IK, format = "(1I0.2)")//SKC_":"// & ! LCOV_EXCL_LINE
                                                    getStr(values(6), length = 2_IK, format = "(1I0.2)")//SKC_":"// & ! LCOV_EXCL_LINE
                                                    getStr(values(7), length = 2_IK, format = "(1I0.2)")//SKC_" "// & ! LCOV_EXCL_LINE
                                                   !getStr(values(7), length = 2_IK, format = "(1I0.2)")//SKC_"."// & ! LCOV_EXCL_LINE
                                                   !getStr(values(8), length = 3_IK, format = "(1I0.3)")//SKC_" "// & ! LCOV_EXCL_LINE
                                                    getStr(values(1))
                elseif (format(i:i) == SKC_"C") then ! Year divided by 100 and truncated to integer (00-99). \warning it works for years up to 9 digits.
                    century = floor(values(1) / 100., IKC)
                    if (abs(century) < 100_IKC) then
                        !write(workspace9(1:3), "(sp,I0.2)") century
                        RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(sp,I0.2)") century
                    else
                        workspace9 = getStr(century)
                        lenSeg = len_trim(workspace9, IKC)
                        RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                        string(epos + 1 : eposnew) = workspace9(1:lenSeg)
                    end if
                elseif (format(i:i) == SKC_"d") then ! Day of the month, zero-padded (01-31).
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") values(3)
                elseif (format(i:i) == SKC_"D") then ! Short MM/DD/YY date, equivalent to %m/%d/%y
                    RESIZE_STRING(8_IK) ! fpp sets eposnew, and resizes string.
                    !write(workspace9(1:2), "(I0.2)") values(2) ! month.
                    !write(workspace9(3:4), "(I0.2)") values(3) ! day.
                    !write(workspace9(5:6), "(I0.2)") mod(abs(values(1)), 100_IKC) ! last two digits of year.
                    !string(epos + 1 : eposnew) =    workspace9(1:2)//sep// & ! LCOV_EXCL_LINE
                    !                                workspace9(3:4)//sep// & ! LCOV_EXCL_LINE
                    !                                workspace9(5:6)
                    write(string(epos + 1 : eposnew), "(I0.2,'/',I0.2,'/',I0.2)") values(2:3), mod(abs(values(1)), 100_IKC) ! last two digits of year.
                elseif (format(i:i) == SKC_"e") then ! Day of the month, zero-padded (01-31).
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I2)") values(3)
                elseif (format(i:i) == SKC_"f") then ! millisecond padded with leading zeros
                    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.3)") values(8)
                elseif (format(i:i) == SKC_"F") then ! Short YYYY-MM-DD date, equivalent to %Y-%m-%d
                    if (values(1) > 0_IKC) then
                        RESIZE_STRING(10_IK) ! fpp sets eposnew, and resizes string.
                    else
                        RESIZE_STRING(11_IK) ! fpp sets eposnew, and resizes string.
                    end if
                    write(string(epos + 1 : eposnew), "(I0.2,'-',I0.2,'-',I0.2)") values(1:3)
                elseif (format(i:i) == SKC_"g") then ! Week-based year, last two digits (00-99).
                    !WeekDate(1:3) = getWeekDate(values(1:3))
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") mod(abs(getWeekYear(values(1:3))), 100_IKC) ! last two digits of the week year.
                elseif (format(i:i) == SKC_"G") then ! Week-based year, full week year, possibly negative.
                    !WeekDate(1:3) = getWeekDate(values(1:3))
                    workspace9 = getStr(getWeekYear(values(1:3))) ! WeekDate(1)) ! full week year, possibly negative.
                    lenSeg = len_trim(workspace9, IKC)
                    RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = workspace9(1:lenSeg)
                elseif (format(i:i) == SKC_"H") then ! Hour in 24h format (00-23)
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") values(5)
                elseif (format(i:i) == SKC_"I") then ! Hour in 12h format (01-12)
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") getHour12(values(5))
                    !if (values(5) < 12_IKC) then
                    !    write(string(epos + 1 : eposnew), "(I0.2)") values(5)
                    !else
                    !    write(string(epos + 1 : eposnew), "(I0.2)") values(5) - 12_IKC
                    !end if
                elseif (format(i:i) == SKC_"j") then ! Day of the year (001-366)
                    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.3)") getOrdinalDay(values(1:3))
                elseif (format(i:i) == SKC_"m") then ! Month as a decimal number (01-12)
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") values(2)
                elseif (format(i:i) == SKC_"M") then ! Minute (00-59)
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") values(6)
                elseif (format(i:i) == SKC_"n") then ! New-line character ('\n')
                    RESIZE_STRING(1_IK) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = achar(10, SKC) ! new_line(SKC_"a")
                elseif (format(i:i) == SKC_"p") then ! AM or PM designation
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    if (values(5) < 12_IKC) then
                        string(epos + 1 : eposnew) = SKC_"AM"
                    else
                        string(epos + 1 : eposnew) = SKC_"PM"
                    end if
                elseif (format(i:i) == SKC_"r") then ! 12-hour clock time *
                    RESIZE_STRING(11_IK) ! fpp sets eposnew, and resizes string.
                    if (values(5) < 12_IKC) then
                        write(string(epos + 1 : eposnew), "(I0.2,':',I0.2,':',I0.2,' am')") getHour12(values(5)), values(6:7)
                    else
                        write(string(epos + 1 : eposnew), "(I0.2,':',I0.2,':',I0.2,' pm')") getHour12(values(5)), values(6:7)
                    end if
                elseif (format(i:i) == SKC_"R") then ! 24-hour HH:MM time, equivalent to %H:%M
                    RESIZE_STRING(5_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2,':',I0.2)") values(5:6)
                elseif (format(i:i) == SKC_"S") then ! Second (00-59)
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") values(7)
                elseif (format(i:i) == SKC_"t") then ! Horizontal-tab character ('\t')
                    RESIZE_STRING(1_IK) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = achar(9, SKC)
                elseif (format(i:i) == SKC_"T") then ! ISO 8601 time format (HH:MM:SS), equivalent to %H:%M:%S
                    RESIZE_STRING(8_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2,':',I0.2,':',I0.2)") values(5:7)
                elseif (format(i:i) == SKC_"u") then ! ISO 8601 weekday as number with Monday as 1 (1-7)
                    RESIZE_STRING(1_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I1)") getWeekDayISO(values(1), values(2), values(3))
                elseif (format(i:i) == SKC_"V") then ! ISO 8601 week number (01-53)
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") getWeekNumber(values(1), values(2), values(3))
                elseif (format(i:i) == SKC_"w") then ! Weekday as a decimal number with Sunday as 0 (0-6)
                    RESIZE_STRING(1_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I1)") getWeekDay(values(1), values(2), values(3))
                elseif (format(i:i) == SKC_"x") then ! Date representation *
                    RESIZE_STRING(8_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2,'/',I0.2,'/',I0.2)") values(2), values(3), mod(abs(values(1)), 100_IKC) ! last two digits of year.
                elseif (format(i:i) == SKC_"X") then ! Time representation *
                    RESIZE_STRING(8_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2,':',I0.2,':',I0.2)") values(5:7)
                elseif (format(i:i) == SKC_"y") then ! Year, last two digits (00-99)
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") mod(abs(values(1)), 100_IKC) ! last two digits of year.
                elseif (format(i:i) == SKC_"Y") then ! Year, last two digits (00-99)
                    workspace9 = getStr(values(1))
                    lenSeg = len_trim(workspace9, IKC)
                    RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = workspace9(1:lenSeg)
                elseif (format(i:i) == SKC_"z") then ! ISO 8601 offset from UTC in timezone in units of minutes
                    RESIZE_STRING(5_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(sp,I0.4)") values(4)
                elseif (format(i:i) == SKC_"Z") then ! Timezone name or abbreviation.
                    abbr = getZoneAbbr(values(4))
                    RESIZE_STRING(len(abbr, IK))
                    string(epos + 1 : eposnew) = abbr
                    !if (values(4) == -60_IKC * 12_IKC) then ! International Day Line West time zone
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"IDLW"
                    !elseif (values(4) == -60_IKC * 11_IKC) then ! Samoa Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"SST"
                    !elseif (values(4) == -60_IKC * 10_IKC) then ! Hawaii–Aleutian Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"HST"
                    !elseif (values(4) == -60_IKC * 9_IKC - 30_IKC) then ! Marquesas Islands Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"MIT"
                    !elseif (values(4) == -60_IKC * 9_IKC) then ! Alaska Standard Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"AKST"
                    !elseif (values(4) == -60_IKC * 8_IKC) then ! Pacific Standard Time (North America)
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"PST"
                    !elseif (values(4) == -60_IKC * 7_IKC) then ! Mountain Standard Time (North America)
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"MST"
                    !elseif (values(4) == -60_IKC * 6_IKC) then ! Central Standard Time (North America)
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"CST"
                    !elseif (values(4) == -60_IKC * 5_IKC) then ! Eastern Standard Time (North America)
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"EST"
                    !elseif (values(4) == -60_IKC * 3_IKC - 30_IKC) then ! Newfoundland Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"NST"
                    !elseif (values(4) == -60_IKC * 3_IKC) then ! Uruguay Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"UYT"
                    !elseif (values(4) == -60_IKC * 2_IKC - 30_IKC) then ! Newfoundland Daylight Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"NDT"
                    !elseif (values(4) == -60_IKC * 2_IKC) then ! Uruguay Summer Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"UYST"
                    !elseif (values(4) == -60_IKC * 1_IKC) then ! Eastern Greenland Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"EGT"
                    !elseif (values(4) == 0_IKC) then ! Coordinated Universal Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"UTC"
                    !elseif (values(4) == +60_IKC * 1_IKC) then ! Central European Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"CET"
                    !elseif (values(4) == +60_IKC * 2_IKC) then ! Eastern European Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"EET"
                    !elseif (values(4) == +60_IKC * 3_IKC) then ! Arabia Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"AST"
                    !elseif (values(4) == +60_IKC * 3_IKC + 30_IKC) then ! Iran Standard Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"IRST"
                    !elseif (values(4) == +60_IKC * 4_IKC) then ! Georgia Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"GET"
                    !elseif (values(4) == +60_IKC * 4_IKC + 30_IKC) then ! Afghanistan Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"AFT"
                    !elseif (values(4) == +60_IKC * 5_IKC) then ! Pakistan Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"PKT"
                    !elseif (values(4) == +60_IKC * 5_IKC + 30_IKC) then ! Indian Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"IST"
                    !elseif (values(4) == +60_IKC * 5_IKC + 45_IKC) then ! Nepal Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"NPT"
                    !elseif (values(4) == +60_IKC * 6_IKC) then ! Bangladesh Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"BST"
                    !elseif (values(4) == +60_IKC * 6_IKC + 30_IKC) then ! Myanmar Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"MMT"
                    !elseif (values(4) == +60_IKC * 7_IKC) then ! Thailand Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"THA"
                    !elseif (values(4) == +60_IKC * 8_IKC) then ! Singapore Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"SST"
                    !elseif (values(4) == +60_IKC * 8_IKC + 45_IKC) then ! Central Western Standard Time (Australia)
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"CWST"
                    !elseif (values(4) == +60_IKC * 9_IKC) then ! Japan Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"JST"
                    !elseif (values(4) == +60_IKC * 9_IKC + 30_IKC) then ! Australian Central Standard Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"ACST"
                    !elseif (values(4) == +60_IKC * 10_IKC) then ! Australian Eastern Standard Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"AEST"
                    !elseif (values(4) == +60_IKC * 10_IKC + 30_IKC) then ! Lord Howe Standard Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"LHST"
                    !elseif (values(4) == +60_IKC * 11_IKC) then ! Pohnpei Standard Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"PONT"
                    !elseif (values(4) == +60_IKC * 12_IKC) then ! New Zealand Standard Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"NZST"
                    !elseif (values(4) == +60_IKC * 12_IKC + 45_IKC) then ! Chatham Standard Time
                    !    RESIZE_STRING(5_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"CHAST"
                    !elseif (values(4) == +60_IKC * 13_IKC) then ! Tonga Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"TOT"
                    !elseif (values(4) == +60_IKC * 13_IKC + 45_IKC) then ! Chatham Daylight Time
                    !    RESIZE_STRING(5_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"CHADT"
                    !elseif (values(4) == +60_IKC * 14_IKC) then ! Line Islands Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKC_"LINT"
                    !end if
                elseif (format(i:i) == SKC_"%") then ! add percentage.
                    RESIZE_STRING(1_IK) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = SKC_"%"
                else ! Unrecognized format.
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = format(i - 1 : i)
                end if
            else ! normal characters.
                RESIZE_STRING(1_IK) ! fpp sets eposnew, and resizes string.
                string(epos + 1 : eposnew) = format(i : i)
            end if
        end do
        if (lenString > eposnew) then
            tempstr = string(1:eposnew)
            call move_alloc(tempstr, string)
        end if
#else
#error  "Unrecognized interface."
#endif

#undef  RESIZE_STRING

        !%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getDateTimeNewZone_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%

        ! The case for no input (current) date and time.
#if     getDateTimeNewZoneC_ENABLED
        integer(IKC) :: values(8)
        call date_and_time(values = values)
        CHECK_ASSERTION(__LINE__, all(values /= -huge(0_IKC)), SK_"@getDateTimeNewZone(): The processor does not have a clock.") ! fpp
        DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, values(1), values(2), values(3), values(4), values(5), values(6), values(7), values(8))
        ! The case for vector of date and time.
#elif   getDateTimeNewZoneV_ENABLED
        integer(IK) :: lenValues
        lenValues = size(values, 1, IK)
        CHECK_ASSERTION(__LINE__, 3_IK < lenValues .and. lenValues < 9_IK, SK_"@getDateTimeNewZone(): The condition `3 < size(values) < 9` must hold. size(values) = "//getStr(size(values))) ! fpp
        if (lenValues == 8_IK) then
            DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, values(1), values(2), values(3), values(4), values(5), values(6), values(7), values(8))
        elseif (lenValues == 7_IK) then
            DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, values(1), values(2), values(3), values(4), values(5), values(6), values(7))
        elseif (lenValues == 6_IK) then
            DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, values(1), values(2), values(3), values(4), values(5), values(6))
        elseif (lenValues == 5_IK) then
            DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, values(1), values(2), values(3), values(4), values(5))
        elseif (lenValues == 4_IK) then
            DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, values(1), values(2), values(3), values(4))
        !elseif (lenValues == 3_IK) then
        !    DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, values(1), values(2), values(3))
        !elseif (lenValues == 2_IK) then
        !    DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, values(1), values(2))
        !elseif (lenValues == 1_IK) then
        !    DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, values(1))
        else
            error stop MODULE_NAME//SK_"@getDateTimeNewZone(): The length of the input argument `values` must be less than 9 and larger than 3. size(values) = "//getStr(lenValues) ! LCOV_EXCL_LINE
        end if
        ! The case for keyword date and time.
#else
#if     !getDateTimeNewZoneYMDZHMSM_ENABLED
        integer(IKC), parameter :: millisecond = 0_IK
#if     !getDateTimeNewZoneYMDZHMS_ENABLED
        integer(IKC), parameter :: second = 0_IK
#if     !getDateTimeNewZoneYMDZHM_ENABLED
        integer(IKC), parameter :: minute = 0_IK
#if     !getDateTimeNewZoneYMDZH_ENABLED
        integer(IKC), parameter :: hour = 0_IK
#if     !getDateTimeNewZoneYMDZ_ENABLED
#error  "Unrecognized interface."
#endif
#endif
#endif
#endif
#endif
        if (zone /= 0_IK .and. newzone /= 0_IK .and. zone /= newzone) then
            DateTimeNewZone(1:8) = getDateTimeUTC(year, month, day, zone, hour, minute, second, millisecond)
            DateTimeNewZone(1:8) = getDateTimeUTC(DateTimeNewZone(1), DateTimeNewZone(2), DateTimeNewZone(3), -newzone, DateTimeNewZone(5), DateTimeNewZone(6), DateTimeNewZone(7), DateTimeNewZone(8))
            DateTimeNewZone(4) = newzone
        elseif (zone /= 0_IK .and. newzone == 0_IK) then
            DateTimeNewZone(1:8) = getDateTimeUTC(year, month, day, zone, hour, minute, second, millisecond)
        elseif (zone == 0_IK .and. newzone /= 0_IK) then
            DateTimeNewZone(1:8) = getDateTimeUTC(year, month, day, -newzone, hour, minute, second, millisecond)
            DateTimeNewZone(4) = newzone
        else
            DateTimeNewZone(1) = year
            DateTimeNewZone(2) = month
            DateTimeNewZone(3) = day
            DateTimeNewZone(4) = zone
            DateTimeNewZone(5) = hour
            DateTimeNewZone(6) = minute
            DateTimeNewZone(7) = second
            DateTimeNewZone(8) = millisecond
        end if
#endif

        !%%%%%%%%%%%%%%%%%%%%%
#elif   getDateTimeUTC_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        ! The case for no input (current) date and time.
#if     getDateTimeUTCC_ENABLED
        integer(IKC) :: values(8)
        call date_and_time(values = values)
        CHECK_ASSERTION(__LINE__, all(values /= -huge(0_IKC)), SK_"@getDateTimeUTC(): The processor does not have a clock.") ! fpp
        DateTimeUTC(1:8) = getDateTimeUTC(values(1), values(2), values(3), values(4), values(5), values(6), values(7), values(8))
        ! The case for vector of date and time.
#elif   getDateTimeUTCV_ENABLED
        integer(IK) :: lenValues
        lenValues = size(values, 1, IK)
        CHECK_ASSERTION(__LINE__, 3_IK < lenValues .and. lenValues < 9_IK, SK_"@getDateTimeUTC(): The condition `3 < size(values) < 9` must hold. size(values) = "//getStr(size(values))) ! fpp
        if (lenValues == 8_IK) then
            DateTimeUTC(1:8) = getDateTimeUTC(values(1), values(2), values(3), values(4), values(5), values(6), values(7), values(8))
        elseif (lenValues == 7_IK) then
            DateTimeUTC(1:8) = getDateTimeUTC(values(1), values(2), values(3), values(4), values(5), values(6), values(7))
        elseif (lenValues == 6_IK) then
            DateTimeUTC(1:8) = getDateTimeUTC(values(1), values(2), values(3), values(4), values(5), values(6))
        elseif (lenValues == 5_IK) then
            DateTimeUTC(1:8) = getDateTimeUTC(values(1), values(2), values(3), values(4), values(5))
        elseif (lenValues == 4_IK) then
            DateTimeUTC(1:8) = getDateTimeUTC(values(1), values(2), values(3), values(4))
        !elseif (lenValues == 3_IK) then
        !    DateTimeUTC(1:8) = getDateTimeUTC(values(1), values(2), values(3))
        !elseif (lenValues == 2_IK) then
        !    DateTimeUTC(1:8) = getDateTimeUTC(values(1), values(2))
        !elseif (lenValues == 1_IK) then
        !    DateTimeUTC(1:8) = getDateTimeUTC(values(1))
        else
            error stop MODULE_NAME//SK_"@getDateTimeUTC(): The length of the input argument `values` must be less than 9 and larger than 3. size(values) = "//getStr(lenValues) ! LCOV_EXCL_LINE
        end if
        ! The case for keyword date and time.
#else
#if     !getDateTimeUTCYMDZHMSM_ENABLED
        integer(IKC), parameter :: millisecond = 0_IK
#if     !getDateTimeUTCYMDZHMS_ENABLED
        integer(IKC), parameter :: second = 0_IK
#if     !getDateTimeUTCYMDZHM_ENABLED
        integer(IKC), parameter :: minute = 0_IK
#if     !getDateTimeUTCYMDZH_ENABLED
        integer(IKC), parameter :: hour = 0_IK
#if     !getDateTimeUTCYMDZ_ENABLED
#error  "Unrecognized interface."
#endif
#endif
#endif
#endif
#endif
        integer(IKC) :: offsetDays
        integer(IKC) :: totalHours
        integer(IKC) :: offsetHours
        integer(IKC) :: totalMinutes
        DateTimeUTC(8) = millisecond
        DateTimeUTC(7) = second
        if (zone /= 0_IKC) then
            totalMinutes = minute - zone
            offsetHours = totalMinutes / 60_IKC
            DateTimeUTC(6) = totalMinutes - offsetHours * 60_IKC
            if (DateTimeUTC(6) < 0_IKC) then
                DateTimeUTC(6) = DateTimeUTC(6) + 60_IKC
                offsetHours = offsetHours - 1_IKC
            end if
            totalHours = hour + offsetHours
            offsetDays = totalHours / 24_IKC
            DateTimeUTC(5) = totalHours - offsetDays * 24_IKC
            if (DateTimeUTC(5) < 0_IKC) then
                DateTimeUTC(5) = DateTimeUTC(5) + 24_IKC
                offsetDays = offsetDays - 1_IKC
            end if
            DateTimeUTC(4) = 0_IKC ! UTC zone
            CHECK_ASSERTION(__LINE__, offsetDays > -2_IKC .and. offsetDays < 2_IKC, SK_"@getDateTimeUTC(): Internal library error: The condition `offsetDays > -2_IKC .and. offsetDays < 2_IKC` must hold. Please report this error to the ParaMonte library developers.")
            if (offsetDays == -1_IKC) then
                DateTimeUTC(1:3) = getDateBefore(year, month, day)
            elseif (offsetDays == +1_IKC) then
                DateTimeUTC(1:3) = getDateAfter(year, month, day)
            else
                DateTimeUTC(3) = day
                DateTimeUTC(2) = month
                DateTimeUTC(1) = year
            end if
        else
            DateTimeUTC(1) = year
            DateTimeUTC(2) = month
            DateTimeUTC(3) = day
            DateTimeUTC(4) = zone
            DateTimeUTC(5) = hour
            DateTimeUTC(6) = minute
            DateTimeUTC(7) = second
            DateTimeUTC(8) = millisecond
        end if
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getDateTimeShifted_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%

        ! The case for no input (current) date and time.
#if     getDateTimeShiftedC_ENABLED
        integer(IKC)    :: values(8)
        call date_and_time(values = values)
        CHECK_ASSERTION(__LINE__, all(values /= -huge(0_IKC)), SK_"@getDateTimeShifted(): The processor does not have a clock.") ! fpp
        dateTimeShifted(1:8) = getDateTimeShifted(amount, values(1), values(2), values(3), values(4), values(5), values(6), values(7), values(8))
        ! The case for vector of date and time.
#elif   getDateTimeShiftedV_ENABLED
        integer(IK) :: lenValues
        lenValues = size(values, 1, IK)
        CHECK_ASSERTION(__LINE__, 0_IK < lenValues .and. lenValues < 9_IK, SK_"@getDateTimeShifted(): The condition `0 < size(values) < 9` must hold. size(values) = "//getStr(size(values))) ! fpp
        if (lenValues == 8_IK) then
            dateTimeShifted(1:8) = getDateTimeShifted(amount, values(1), values(2), values(3), values(4), values(5), values(6), values(7), values(8))
        elseif (lenValues == 7_IK) then
            dateTimeShifted(1:8) = getDateTimeShifted(amount, values(1), values(2), values(3), values(4), values(5), values(6), values(7))
        elseif (lenValues == 6_IK) then
            dateTimeShifted(1:8) = getDateTimeShifted(amount, values(1), values(2), values(3), values(4), values(5), values(6))
        elseif (lenValues == 5_IK) then
            dateTimeShifted(1:8) = getDateTimeShifted(amount, values(1), values(2), values(3), values(4), values(5))
        elseif (lenValues == 4_IK) then
            dateTimeShifted(1:8) = getDateTimeShifted(amount, values(1), values(2), values(3), values(4))
        elseif (lenValues == 3_IK) then
            dateTimeShifted(1:8) = getDateTimeShifted(amount, values(1), values(2), values(3))
        elseif (lenValues == 2_IK) then
            dateTimeShifted(1:8) = getDateTimeShifted(amount, values(1), values(2))
        elseif (lenValues == 1_IK) then
            dateTimeShifted(1:8) = getDateTimeShifted(amount, values(1))
        else
            error stop MODULE_NAME//SK_"@getDateTimeShifted(): The length of the input argument `values` must be less than 9 and non-zero. size(values) = "//getStr(lenValues) ! LCOV_EXCL_LINE
        end if
        ! The case for keyword date and time.
#else
#if     !getDateTimeShiftedYMDZHMSM_ENABLED
        integer(IKC), parameter :: millisecond = 0_IK
#if     !getDateTimeShiftedYMDZHMS_ENABLED
        integer(IKC), parameter :: second = 0_IK
#if     !getDateTimeShiftedYMDZHM_ENABLED
        integer(IKC), parameter :: minute = 0_IK
#if     !getDateTimeShiftedYMDZH_ENABLED
        integer(IKC), parameter :: hour = 0_IK
#if     !getDateTimeShiftedYMDZ_ENABLED
        integer(IKC), parameter :: zone = 0_IK
#if     !getDateTimeShiftedYMD_ENABLED
        integer(IKC), parameter :: day = 1_IK
#if     !getDateTimeShiftedYM_ENABLED
        integer(IKC), parameter :: month = 1_IK
#if     !getDateTimeShiftedY_ENABLED
#error  "Unrecognized interface."
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
        real(RKC) :: julianDay
        CHECK_ASSERTION(__LINE__, isValidDateTime(year, month, day, zone, hour, minute, second, millisecond), \
        SK_"@getDateTimeShifted(): The specified Gregorian date and time must be valid and consistent. [year, month, day, zone, hour, minute, second, millisecond] = "// \
        getStr([year, month, day, zone, hour, minute, second, millisecond])) ! fpp
        julianDay = getJulianDay(year, month, day, 0_IK, hour, minute, second, millisecond)
        !dateTimeShifted(1:8) = getDateTime(julianDay + amount, zone)
        dateTimeShifted(1:8) = getDateTime(julianDay + amount)
        dateTimeShifted(4) = zone
#endif

        !%%%%%%%%%%%%%%%%%%%
#elif   getJulianDay_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        ! The case for no input (current) date and time.
#if     getJulianDayC_ENABLED
        integer(IKC) :: values(8)
        call date_and_time(values = values)
        CHECK_ASSERTION(__LINE__, all(values /= -huge(0_IKC)), SK_"@getJulianDay(): The processor does not have a clock.") ! fpp
        julianDay = getJulianDay(values(1), values(2), values(3), values(4), values(5), values(6), values(7), values(8))
        ! The case for vector of date and time.
#elif   getJulianDayV_ENABLED
        integer(IK) :: lenValues
        lenValues = size(values, 1, IK)
        CHECK_ASSERTION(__LINE__, 0_IK < lenValues .and. lenValues < 9_IK, SK_"@getJulianDay(): The condition `0 < size(values) < 9` must hold. size(values) = "//getStr(size(values))) ! fpp
        if (lenValues == 8_IK) then
            julianDay = getJulianDay(values(1), values(2), values(3), values(4), values(5), values(6), values(7), values(8))
        elseif (lenValues == 7_IK) then
            julianDay = getJulianDay(values(1), values(2), values(3), values(4), values(5), values(6), values(7))
        elseif (lenValues == 6_IK) then
            julianDay = getJulianDay(values(1), values(2), values(3), values(4), values(5), values(6))
        elseif (lenValues == 5_IK) then
            julianDay = getJulianDay(values(1), values(2), values(3), values(4), values(5))
        elseif (lenValues == 4_IK) then
            julianDay = getJulianDay(values(1), values(2), values(3), values(4))
        elseif (lenValues == 3_IK) then
            julianDay = getJulianDay(values(1), values(2), values(3))
        elseif (lenValues == 2_IK) then
            julianDay = getJulianDay(values(1), values(2))
        elseif (lenValues == 1_IK) then
            julianDay = getJulianDay(values(1))
        else
            error stop MODULE_NAME//SK_"@getJulianDay(): The length of the input argument `values` must be less than 9 and non-zero. size(values) = "//getStr(lenValues) ! LCOV_EXCL_LINE
        end if
        ! The case for keyword date and time.
#else
#if     !getJulianDayYMDZHMSM_ENABLED
        integer(IKC), parameter :: millisecond = 0_IK
#if     !getJulianDayYMDZHMS_ENABLED
        integer(IKC), parameter :: second = 0_IK
#if     !getJulianDayYMDZHM_ENABLED
        integer(IKC), parameter :: minute = 0_IK
#if     !getJulianDayYMDZH_ENABLED
        integer(IKC), parameter :: hour = 0_IK
#if     !getJulianDayYMDZ_ENABLED
        integer(IKC), parameter :: zone = 0_IK
#if     !getJulianDayYMD_ENABLED
        integer(IKC), parameter :: day = 1_IK
#if     !getJulianDayYM_ENABLED
        integer(IKC), parameter :: month = 1_IK
#if     !getJulianDayY_ENABLED
#error  "Unrecognized interface."
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
        real(RKC)   , parameter :: DAY_PER_HOUR = 1._RKC / 24._RKC
        real(RKC)   , parameter :: DAY_PER_MIN = 1._RKC / 1440._RKC
        real(RKC)   , parameter :: DAY_PER_SEC = 1._RKC / 86400._RKC
        integer(IKC), parameter :: VECTOR(12) = [306_IKC, 337_IKC, 0_IKC, 31_IKC, 61_IKC, 92_IKC, 122_IKC, 153_IKC, 184_IKC, 214_IKC, 245_IKC, 275_IKC]
        integer(IKC) :: yearCorrected, DateTimeUTC(8)

        CHECK_ASSERTION(__LINE__, isValidDateTime(year, month, day, zone, hour, minute, second, millisecond), \
        SK_"@getJulianDay(): The specified Gregorian date and time must be valid and consistent. [year, month, day, zone, hour, minute, second, millisecond] = "// \
        getStr([year, month, day, zone, hour, minute, second, millisecond])) ! fpp

        if (zone /= 0_IK) then
            DateTimeUTC(1:8) = getDateTimeUTC(year, month, day, zone, hour, minute, second, millisecond)
        else
            DateTimeUTC(1) = year
            DateTimeUTC(2) = month
            DateTimeUTC(3) = day
            DateTimeUTC(4) = zone
            DateTimeUTC(5) = hour
            DateTimeUTC(6) = minute
            DateTimeUTC(7) = second
            DateTimeUTC(8) = millisecond
        end if

        if (DateTimeUTC(2) == 1_IKC .or. DateTimeUTC(2) == 2_IKC) then
            yearCorrected = DateTimeUTC(1) - 1_IKC
        else
            yearCorrected = DateTimeUTC(1)
        end if
        julianDay   = 1721118.5_RKC & ! LCOV_EXCL_LINE
                    + real(DateTimeUTC(3), RKC) & ! LCOV_EXCL_LINE
                    + real(VECTOR(DateTimeUTC(2)), RKC) & ! LCOV_EXCL_LINE
                    + 365._RKC * real(yearCorrected, RKC) & ! LCOV_EXCL_LINE
                    + real(floor(yearCorrected * 0.25_RKC, IKC) - floor(yearCorrected * 0.01_RKC, IKC) + floor(yearCorrected * 0.0025_RKC, IKC), RKC) & ! LCOV_EXCL_LINE
                    + DateTimeUTC(5) * DAY_PER_HOUR + DateTimeUTC(6) * DAY_PER_MIN + (real(DateTimeUTC(7),RKC) + 0.001_RKC * DateTimeUTC(8)) * DAY_PER_SEC
#endif

        !%%%%%%%%%%%%%%%%%%%%%%
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%
