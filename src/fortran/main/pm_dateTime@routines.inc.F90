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
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, March 23, 2012, 5:33 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%
#if     isValidZone_ENABLED
        !%%%%%%%%%%%%%%%%%%

        isValid = logical(int(ZONE_MIN, IKG) <= zone .and. zone <= int(ZONE_MAX, IKG), LK)

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
        isValid = .true._LK ! year /= 0_IKG ! year zero is explicitly allowed in ISO 8601.
#elif   isValidDateTimeYM_ENABLED
        isValid = isValidDateTime(year) .and. 0_IKG < month .and. month < 13_IKG
#elif   isValidDateTimeYMD_ENABLED
        isValid = isValidDateTime(year, month) .and. day > 0_IKG
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
        isValid = isValidDateTime(year, month, day, zone) .and. 0_IKG <= hour .and. hour < 24_IKG
#elif   isValidDateTimeYMDZHM_ENABLED
        isValid = isValidDateTime(year, month, day, zone, hour) .and. 0_IKG <= minute .and. minute < 60_IKG
#elif   isValidDateTimeYMDZHMS_ENABLED
        isValid = isValidDateTime(year, month, day, zone, hour, minute) .and. 0_IKG <= second .and. second < 60_IKG
#elif   isValidDateTimeYMDZHMSM_ENABLED
        isValid = isValidDateTime(year, month, day, zone, hour, minute, second) .and. 0_IKG <= millisecond .and. millisecond < 1000_IKG
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getDateTimeValues_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

#if     getDateTimeValuesJ_ENABLED
        real(RKG)   , parameter :: K1 = 0.25_RKG ! The condition `0.002929687499688476 < K1 <= 0.2521972656249999` must hold.
        real(RKG)   , parameter :: K2 = 0.25_RKG ! The condition `0. < K2 <= 0.25` must hold.
        real(RKG)   , parameter :: MEAN_YEAR = 365.25_RKG
        real(RKG)   , parameter :: MEAN_YEAR_INVERSE = 1._RKG / MEAN_YEAR
        real(RKG)   , parameter :: WHOLE_CENTURY_FACTOR = 1._RKG / 36524.25_RKG
        integer(IKG), parameter :: DAYS_IN_PAST_MONTHS(3:14) = [0_IKG, 31_IKG, 61_IKG, 92_IKG, 122_IKG, 153_IKG, 184_IKG, 214_IKG, 245_IKG, 275_IKG, 306_IKG, 337_IKG] ! last two correspond to jan, feb.
        integer(IKG)    :: julianDayOffset_IK ! Z
        integer(IKG)    :: fullCenturyCount
        integer(IKG)    :: dayOfYearTerm
        real(RKG)       :: hours_RK
        real(RKG)       :: minutes_RK
        real(RKG)       :: seconds_RK
        real(RKG)       :: julianDayOffset_RK
        real(RKG)       :: julianDayOffsetResidual ! R
        real(RKG)       :: daysInWholeCenturyMinusConst ! B
        julianDayOffset_RK = julianDay - 1721118.5_RKG
        julianDayOffset_IK = floor(julianDayOffset_RK, IKG) ! Z
        julianDayOffsetResidual = julianDayOffset_RK - real(julianDayOffset_IK, RKG) ! R : always positive
        fullCenturyCount = floor((real(julianDayOffset_IK, RKG) - K1) * WHOLE_CENTURY_FACTOR, IKG) ! A
        dayOfYearTerm = julianDayOffset_IK + fullCenturyCount - floor(fullCenturyCount * 0.25_RKG, IKG) ! Z + A - floor(A/4)
        daysInWholeCenturyMinusConst = dayOfYearTerm - K2 ! B
        values(1) = floor(daysInWholeCenturyMinusConst * MEAN_YEAR_INVERSE, IKG) ! Y : Calendar Year Starting March.
        values(3) = dayOfYearTerm - floor(values(1) * MEAN_YEAR, IKG) ! C : Day Of Year.
        values(2) = (5_IKG * values(3) + 456_IKG) / 153_IKG ! M : Month Of Year in the range 3:14.
        values(3) = values(3) - DAYS_IN_PAST_MONTHS(values(2)) !+ int(julianDayOffsetResidual, IKG)
        !if (julianDayOffsetResidual > 1._RKG) error stop getStr(julianDayOffsetResidual)
        if (values(2) > 12_IKG) then
            values(1) = values(1) + 1_IKG
            values(2) = values(2) - 12_IKG
        end if
        values(4) = 0_IKG
        hours_RK = julianDayOffsetResidual * 24._RKG
        values(5) = int(hours_RK, IKG)
        !if (values(5) < 0_IKG) then
        !    values(1:3) = getDateBefore(values(1), values(2), values(3))
        !    values(5) = 24_IKG - values(5)
        !end if
        minutes_RK = (hours_RK - values(5)) * 60._RKG
        values(6) = int(minutes_RK, IKG)
        seconds_RK = (minutes_RK - values(6)) * 60._RKG
        values(7) = int(seconds_RK, IKG)
        values(8) = nint((seconds_RK - values(7)) * 1000._RKG, IKG)
        if (values(8) == 1000_IKG) then
            values(8) = 0_IKG
            values(7) = values(7) + 1_IKG
            if (values(7) == 60_IKG) then
                values(7) = 0_IKG
                values(6) = values(6) + 1_IKG
                if (values(6) == 60_IKG) then
                    values(6) = 0_IKG
                    values(5) = values(5) + 1_IKG
                end if
                if (values(5) == 24_IKG) then
                    values(5) = 0_IKG
                    values(1:3) = getDateAfter(values(1), values(2), values(3))
                end if
            end if
        end if
#elif   getDateTimeValuesJZ_ENABLED
        values(1:8) = getDateTimeNewZone(zone, getDateTime(julianDay))
#elif   getDateTimeValuesC_ENABLED
        call date_and_time(values = values)
        !integer(IK)     :: lenValues
        !integer(IKG)    :: Values_(8)
        !lenValues = size(values, 1, IK)
        !check_assertion(__LINE__, 0_IK < lenValues .and. lenValues < 9_IK, SK_"@getDateTime(): The input argument `values` must have a non-zero size that is less than 9. size(values) = "//getStr(lenValues)) ! fpp
        !if (lenValues > 7_IK) then
        !    call date_and_time(values = values)
        !else
        !    call date_and_time(values = Values_)
        !    values(1:lenValues) = Values_(1:lenValues)
        !end if
        !check_assertion(__LINE__, all(values(1:lenValues) /= -huge(0_IKG)), SK_"@getDateTime(): The processor does not have a clock.") ! fpp
#elif   getDateTimeValuesY_ENABLED
        values(1) = year
        values(2) = 1_IKG
        values(3) = 1_IKG
        values(4:8) = 0_IKG
#elif   getDateTimeValuesYM_ENABLED
        values(1) = year
        values(2) = month
        values(3) = 1_IKG
        values(4:8) = 0_IKG
#elif   getDateTimeValuesYMD_ENABLED
        values(1) = year
        values(2) = month
        values(3) = day
        values(4:8) = 0_IKG
#elif   getDateTimeValuesYMDZ_ENABLED
        values(1) = year
        values(2) = month
        values(3) = day
        values(4) = zone
        values(5:8) = 0_IKG
#elif   getDateTimeValuesYMDZH_ENABLED
        values(1) = year
        values(2) = month
        values(3) = day
        values(4) = zone
        values(5) = hour
        values(6:8) = 0_IKG
#elif   getDateTimeValuesYMDZHM_ENABLED
        values(1) = year
        values(2) = month
        values(3) = day
        values(4) = zone
        values(5) = hour
        values(6) = minute
        values(7:8) = 0_IKG
#elif   getDateTimeValuesYMDZHMS_ENABLED
        values(1) = year
        values(2) = month
        values(3) = day
        values(4) = zone
        values(5) = hour
        values(6) = minute
        values(7) = second
        values(8) = 0_IKG
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
allocate(character(eposnew,SKG) :: tempstr); \
tempstr(1:lenString) = string; \
call move_alloc(tempstr, string); \
lenString = eposnew; \
end if;

#if     getDateTimeStringC_ENABLED
        integer(IKG)    :: values(8)
        call date_and_time(values = values)
        string = getDateTime(format, values)
#elif   getDateTimeStringV_ENABLED
        !>  \warning
        !>  The output of getStr() in this procedure is of kind \SK which is incompatible with any value of SKG /= SK.
        !>  For now, this is not an issue since both kinds point to the default character kind.
        !>  This will however become an issue once Fortran standard and compilers support non-default date and time characters.
        use pm_val2str, only: getStr
        character(:,SKG), allocatable   :: tempstr, abbr
       !character(1,SKG)                :: sep
       !character(28,SKG)               :: workspace
        character(9,SKG)                :: workspace9
        integer(IK)                     :: lenString, lenFormat, i, epos, eposnew, lenSeg
        integer(IKG)                    :: century!, WeekDate(3)
        allocate(character(127,SKG)     :: string)
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
            if (format(i:i) == SKG_"%") then
                if (i == lenFormat) exit
                i = i + 1_IK
                if (format(i:i) == SKG_"a") then ! Abbreviated weekday name *
                    lenSeg = 3_IK
                    RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = WEEKDAY_NAME_ISO(getWeekDayISO(values(1), values(2), values(3)))(1:lenSeg)
                elseif (format(i:i) == SKG_"A") then ! Full weekday name *
                    workspace9 = WEEKDAY_NAME_ISO(getWeekDayISO(values(1), values(2), values(3)))
                    lenSeg = len_trim(workspace9, IKG)
                    RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = workspace9(1:lenSeg)
                elseif (format(i:i) == SKG_"b" .or. format(i:i) == SKG_"h") then ! Abbreviated month name * .or. Abbreviated month name * (same as %b)
                    lenSeg = 3_IK
                    RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = MONTH_NAME(values(2))(1:lenSeg)
                elseif (format(i:i) == SKG_"B") then ! Full month name *
                    lenSeg = len_trim(MONTH_NAME(values(2)), IKG)
                    RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = MONTH_NAME(values(2))(1:lenSeg)
                elseif (format(i:i) == SKG_"c") then ! Date and time representation *
                    if (values(1) > 0_IKG) then
                        RESIZE_STRING(24_IK) ! fpp sets eposnew, and resizes string.
                    else
                        RESIZE_STRING(25_IK) ! fpp sets eposnew, and resizes string.
                    end if
                    string(epos + 1 : eposnew) =    WEEKDAY_NAME_ISO(getWeekDayISO(values(1), values(2), values(3)))(1:3)//SKG_" "// & ! LCOV_EXCL_LINE
                                                    MONTH_NAME(values(2))(1:3)//SKG_" "// & ! LCOV_EXCL_LINE
                                                    getStr(values(3), length = 2_IK, format = "(1I0.2)")//SKG_" "// & ! LCOV_EXCL_LINE
                                                    getStr(values(5), length = 2_IK, format = "(1I0.2)")//SKG_":"// & ! LCOV_EXCL_LINE
                                                    getStr(values(6), length = 2_IK, format = "(1I0.2)")//SKG_":"// & ! LCOV_EXCL_LINE
                                                    getStr(values(7), length = 2_IK, format = "(1I0.2)")//SKG_" "// & ! LCOV_EXCL_LINE
                                                   !getStr(values(7), length = 2_IK, format = "(1I0.2)")//SKG_"."// & ! LCOV_EXCL_LINE
                                                   !getStr(values(8), length = 3_IK, format = "(1I0.3)")//SKG_" "// & ! LCOV_EXCL_LINE
                                                    getStr(values(1))
                elseif (format(i:i) == SKG_"C") then ! Year divided by 100 and truncated to integer (00-99). \warning it works for years up to 9 digits.
                    century = floor(values(1) / 100., IKG)
                    if (abs(century) < 100_IKG) then
                        !write(workspace9(1:3), "(sp,I0.2)") century
                        RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(sp,I0.2)") century
                    else
                        workspace9 = getStr(century)
                        lenSeg = len_trim(workspace9, IKG)
                        RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                        string(epos + 1 : eposnew) = workspace9(1:lenSeg)
                    end if
                elseif (format(i:i) == SKG_"d") then ! Day of the month, zero-padded (01-31).
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") values(3)
                elseif (format(i:i) == SKG_"D") then ! Short MM/DD/YY date, equivalent to %m/%d/%y
                    RESIZE_STRING(8_IK) ! fpp sets eposnew, and resizes string.
                    !write(workspace9(1:2), "(I0.2)") values(2) ! month.
                    !write(workspace9(3:4), "(I0.2)") values(3) ! day.
                    !write(workspace9(5:6), "(I0.2)") mod(abs(values(1)), 100_IKG) ! last two digits of year.
                    !string(epos + 1 : eposnew) =    workspace9(1:2)//sep// & ! LCOV_EXCL_LINE
                    !                                workspace9(3:4)//sep// & ! LCOV_EXCL_LINE
                    !                                workspace9(5:6)
                    write(string(epos + 1 : eposnew), "(I0.2,'/',I0.2,'/',I0.2)") values(2:3), mod(abs(values(1)), 100_IKG) ! last two digits of year.
                elseif (format(i:i) == SKG_"e") then ! Day of the month, zero-padded (01-31).
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I2)") values(3)
                elseif (format(i:i) == SKG_"f") then ! millisecond padded with leading zeros
                    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.3)") values(8)
                elseif (format(i:i) == SKG_"F") then ! Short YYYY-MM-DD date, equivalent to %Y-%m-%d
                    if (values(1) > 0_IKG) then
                        RESIZE_STRING(10_IK) ! fpp sets eposnew, and resizes string.
                    else
                        RESIZE_STRING(11_IK) ! fpp sets eposnew, and resizes string.
                    end if
                    write(string(epos + 1 : eposnew), "(I0.2,'-',I0.2,'-',I0.2)") values(1:3)
                elseif (format(i:i) == SKG_"g") then ! Week-based year, last two digits (00-99).
                    !WeekDate(1:3) = getWeekDate(values(1:3))
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") mod(abs(getWeekYear(values(1:3))), 100_IKG) ! last two digits of the week year.
                elseif (format(i:i) == SKG_"G") then ! Week-based year, full week year, possibly negative.
                    !WeekDate(1:3) = getWeekDate(values(1:3))
                    workspace9 = getStr(getWeekYear(values(1:3))) ! WeekDate(1)) ! full week year, possibly negative.
                    lenSeg = len_trim(workspace9, IKG)
                    RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = workspace9(1:lenSeg)
                elseif (format(i:i) == SKG_"H") then ! Hour in 24h format (00-23)
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") values(5)
                elseif (format(i:i) == SKG_"I") then ! Hour in 12h format (01-12)
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") getHour12(values(5))
                    !if (values(5) < 12_IKG) then
                    !    write(string(epos + 1 : eposnew), "(I0.2)") values(5)
                    !else
                    !    write(string(epos + 1 : eposnew), "(I0.2)") values(5) - 12_IKG
                    !end if
                elseif (format(i:i) == SKG_"j") then ! Day of the year (001-366)
                    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.3)") getOrdinalDay(values(1:3))
                elseif (format(i:i) == SKG_"m") then ! Month as a decimal number (01-12)
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") values(2)
                elseif (format(i:i) == SKG_"M") then ! Minute (00-59)
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") values(6)
                elseif (format(i:i) == SKG_"n") then ! New-line character ('\n')
                    RESIZE_STRING(1_IK) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = achar(10, SKG) ! new_line(SKG_"a")
                elseif (format(i:i) == SKG_"p") then ! AM or PM designation
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    if (values(5) < 12_IKG) then
                        string(epos + 1 : eposnew) = SKG_"AM"
                    else
                        string(epos + 1 : eposnew) = SKG_"PM"
                    end if
                elseif (format(i:i) == SKG_"r") then ! 12-hour clock time *
                    RESIZE_STRING(11_IK) ! fpp sets eposnew, and resizes string.
                    if (values(5) < 12_IKG) then
                        write(string(epos + 1 : eposnew), "(I0.2,':',I0.2,':',I0.2,' am')") getHour12(values(5)), values(6:7)
                    else
                        write(string(epos + 1 : eposnew), "(I0.2,':',I0.2,':',I0.2,' pm')") getHour12(values(5)), values(6:7)
                    end if
                elseif (format(i:i) == SKG_"R") then ! 24-hour HH:MM time, equivalent to %H:%M
                    RESIZE_STRING(5_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2,':',I0.2)") values(5:6)
                elseif (format(i:i) == SKG_"S") then ! Second (00-59)
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") values(7)
                elseif (format(i:i) == SKG_"t") then ! Horizontal-tab character ('\t')
                    RESIZE_STRING(1_IK) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = achar(9, SKG)
                elseif (format(i:i) == SKG_"T") then ! ISO 8601 time format (HH:MM:SS), equivalent to %H:%M:%S
                    RESIZE_STRING(8_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2,':',I0.2,':',I0.2)") values(5:7)
                elseif (format(i:i) == SKG_"u") then ! ISO 8601 weekday as number with Monday as 1 (1-7)
                    RESIZE_STRING(1_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I1)") getWeekDayISO(values(1), values(2), values(3))
                elseif (format(i:i) == SKG_"V") then ! ISO 8601 week number (01-53)
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") getWeekNumber(values(1), values(2), values(3))
                elseif (format(i:i) == SKG_"w") then ! Weekday as a decimal number with Sunday as 0 (0-6)
                    RESIZE_STRING(1_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I1)") getWeekDay(values(1), values(2), values(3))
                elseif (format(i:i) == SKG_"x") then ! Date representation *
                    RESIZE_STRING(8_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2,'/',I0.2,'/',I0.2)") values(2), values(3), mod(abs(values(1)), 100_IKG) ! last two digits of year.
                elseif (format(i:i) == SKG_"X") then ! Time representation *
                    RESIZE_STRING(8_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2,':',I0.2,':',I0.2)") values(5:7)
                elseif (format(i:i) == SKG_"y") then ! Year, last two digits (00-99)
                    RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(I0.2)") mod(abs(values(1)), 100_IKG) ! last two digits of year.
                elseif (format(i:i) == SKG_"Y") then ! Year, last two digits (00-99)
                    workspace9 = getStr(values(1))
                    lenSeg = len_trim(workspace9, IKG)
                    RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = workspace9(1:lenSeg)
                elseif (format(i:i) == SKG_"z") then ! ISO 8601 offset from UTC in timezone in units of minutes
                    RESIZE_STRING(5_IK) ! fpp sets eposnew, and resizes string.
                    write(string(epos + 1 : eposnew), "(sp,I0.4)") values(4)
                elseif (format(i:i) == SKG_"Z") then ! Timezone name or abbreviation.
                    abbr = getZoneAbbr(values(4))
                    RESIZE_STRING(len(abbr, IK))
                    string(epos + 1 : eposnew) = abbr
                    !if (values(4) == -60_IKG * 12_IKG) then ! International Day Line West time zone
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"IDLW"
                    !elseif (values(4) == -60_IKG * 11_IKG) then ! Samoa Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"SST"
                    !elseif (values(4) == -60_IKG * 10_IKG) then ! Hawaii–Aleutian Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"HST"
                    !elseif (values(4) == -60_IKG * 9_IKG - 30_IKG) then ! Marquesas Islands Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"MIT"
                    !elseif (values(4) == -60_IKG * 9_IKG) then ! Alaska Standard Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"AKST"
                    !elseif (values(4) == -60_IKG * 8_IKG) then ! Pacific Standard Time (North America)
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"PST"
                    !elseif (values(4) == -60_IKG * 7_IKG) then ! Mountain Standard Time (North America)
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"MST"
                    !elseif (values(4) == -60_IKG * 6_IKG) then ! Central Standard Time (North America)
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"CST"
                    !elseif (values(4) == -60_IKG * 5_IKG) then ! Eastern Standard Time (North America)
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"EST"
                    !elseif (values(4) == -60_IKG * 3_IKG - 30_IKG) then ! Newfoundland Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"NST"
                    !elseif (values(4) == -60_IKG * 3_IKG) then ! Uruguay Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"UYT"
                    !elseif (values(4) == -60_IKG * 2_IKG - 30_IKG) then ! Newfoundland Daylight Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"NDT"
                    !elseif (values(4) == -60_IKG * 2_IKG) then ! Uruguay Summer Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"UYST"
                    !elseif (values(4) == -60_IKG * 1_IKG) then ! Eastern Greenland Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"EGT"
                    !elseif (values(4) == 0_IKG) then ! Coordinated Universal Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"UTC"
                    !elseif (values(4) == +60_IKG * 1_IKG) then ! Central European Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"CET"
                    !elseif (values(4) == +60_IKG * 2_IKG) then ! Eastern European Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"EET"
                    !elseif (values(4) == +60_IKG * 3_IKG) then ! Arabia Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"AST"
                    !elseif (values(4) == +60_IKG * 3_IKG + 30_IKG) then ! Iran Standard Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"IRST"
                    !elseif (values(4) == +60_IKG * 4_IKG) then ! Georgia Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"GET"
                    !elseif (values(4) == +60_IKG * 4_IKG + 30_IKG) then ! Afghanistan Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"AFT"
                    !elseif (values(4) == +60_IKG * 5_IKG) then ! Pakistan Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"PKT"
                    !elseif (values(4) == +60_IKG * 5_IKG + 30_IKG) then ! Indian Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"IST"
                    !elseif (values(4) == +60_IKG * 5_IKG + 45_IKG) then ! Nepal Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"NPT"
                    !elseif (values(4) == +60_IKG * 6_IKG) then ! Bangladesh Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"BST"
                    !elseif (values(4) == +60_IKG * 6_IKG + 30_IKG) then ! Myanmar Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"MMT"
                    !elseif (values(4) == +60_IKG * 7_IKG) then ! Thailand Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"THA"
                    !elseif (values(4) == +60_IKG * 8_IKG) then ! Singapore Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"SST"
                    !elseif (values(4) == +60_IKG * 8_IKG + 45_IKG) then ! Central Western Standard Time (Australia)
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"CWST"
                    !elseif (values(4) == +60_IKG * 9_IKG) then ! Japan Standard Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"JST"
                    !elseif (values(4) == +60_IKG * 9_IKG + 30_IKG) then ! Australian Central Standard Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"ACST"
                    !elseif (values(4) == +60_IKG * 10_IKG) then ! Australian Eastern Standard Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"AEST"
                    !elseif (values(4) == +60_IKG * 10_IKG + 30_IKG) then ! Lord Howe Standard Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"LHST"
                    !elseif (values(4) == +60_IKG * 11_IKG) then ! Pohnpei Standard Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"PONT"
                    !elseif (values(4) == +60_IKG * 12_IKG) then ! New Zealand Standard Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"NZST"
                    !elseif (values(4) == +60_IKG * 12_IKG + 45_IKG) then ! Chatham Standard Time
                    !    RESIZE_STRING(5_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"CHAST"
                    !elseif (values(4) == +60_IKG * 13_IKG) then ! Tonga Time
                    !    RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"TOT"
                    !elseif (values(4) == +60_IKG * 13_IKG + 45_IKG) then ! Chatham Daylight Time
                    !    RESIZE_STRING(5_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"CHADT"
                    !elseif (values(4) == +60_IKG * 14_IKG) then ! Line Islands Time
                    !    RESIZE_STRING(4_IK) ! fpp sets eposnew, and resizes string.
                    !    string(epos + 1 : eposnew) = SKG_"LINT"
                    !end if
                elseif (format(i:i) == SKG_"%") then ! add percentage.
                    RESIZE_STRING(1_IK) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = SKG_"%"
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
        integer(IKG) :: values(8)
        call date_and_time(values = values)
        CHECK_ASSERTION(__LINE__, all(values /= -huge(0_IKG)), SK_"@getDateTimeNewZone(): The processor does not have a clock.") ! fpp
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
        integer(IKG), parameter :: millisecond = 0_IK
#if     !getDateTimeNewZoneYMDZHMS_ENABLED
        integer(IKG), parameter :: second = 0_IK
#if     !getDateTimeNewZoneYMDZHM_ENABLED
        integer(IKG), parameter :: minute = 0_IK
#if     !getDateTimeNewZoneYMDZH_ENABLED
        integer(IKG), parameter :: hour = 0_IK
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
        integer(IKG) :: values(8)
        call date_and_time(values = values)
        CHECK_ASSERTION(__LINE__, all(values /= -huge(0_IKG)), SK_"@getDateTimeUTC(): The processor does not have a clock.") ! fpp
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
        integer(IKG), parameter :: millisecond = 0_IK
#if     !getDateTimeUTCYMDZHMS_ENABLED
        integer(IKG), parameter :: second = 0_IK
#if     !getDateTimeUTCYMDZHM_ENABLED
        integer(IKG), parameter :: minute = 0_IK
#if     !getDateTimeUTCYMDZH_ENABLED
        integer(IKG), parameter :: hour = 0_IK
#if     !getDateTimeUTCYMDZ_ENABLED
#error  "Unrecognized interface."
#endif
#endif
#endif
#endif
#endif
        integer(IKG) :: offsetDays
        integer(IKG) :: totalHours
        integer(IKG) :: offsetHours
        integer(IKG) :: totalMinutes
        DateTimeUTC(8) = millisecond
        DateTimeUTC(7) = second
        if (zone /= 0_IKG) then
            totalMinutes = minute - zone
            offsetHours = totalMinutes / 60_IKG
            DateTimeUTC(6) = totalMinutes - offsetHours * 60_IKG
            if (DateTimeUTC(6) < 0_IKG) then
                DateTimeUTC(6) = DateTimeUTC(6) + 60_IKG
                offsetHours = offsetHours - 1_IKG
            end if
            totalHours = hour + offsetHours
            offsetDays = totalHours / 24_IKG
            DateTimeUTC(5) = totalHours - offsetDays * 24_IKG
            if (DateTimeUTC(5) < 0_IKG) then
                DateTimeUTC(5) = DateTimeUTC(5) + 24_IKG
                offsetDays = offsetDays - 1_IKG
            end if
            DateTimeUTC(4) = 0_IKG ! UTC zone
            CHECK_ASSERTION(__LINE__, offsetDays > -2_IKG .and. offsetDays < 2_IKG, SK_"@getDateTimeUTC(): Internal library error: The condition `offsetDays > -2_IKG .and. offsetDays < 2_IKG` must hold. Please report this error to the ParaMonte library developers.")
            if (offsetDays == -1_IKG) then
                DateTimeUTC(1:3) = getDateBefore(year, month, day)
            elseif (offsetDays == +1_IKG) then
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
        integer(IKG)    :: values(8)
        call date_and_time(values = values)
        CHECK_ASSERTION(__LINE__, all(values /= -huge(0_IKG)), SK_"@getDateTimeShifted(): The processor does not have a clock.") ! fpp
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
        integer(IKG), parameter :: millisecond = 0_IK
#if     !getDateTimeShiftedYMDZHMS_ENABLED
        integer(IKG), parameter :: second = 0_IK
#if     !getDateTimeShiftedYMDZHM_ENABLED
        integer(IKG), parameter :: minute = 0_IK
#if     !getDateTimeShiftedYMDZH_ENABLED
        integer(IKG), parameter :: hour = 0_IK
#if     !getDateTimeShiftedYMDZ_ENABLED
        integer(IKG), parameter :: zone = 0_IK
#if     !getDateTimeShiftedYMD_ENABLED
        integer(IKG), parameter :: day = 1_IK
#if     !getDateTimeShiftedYM_ENABLED
        integer(IKG), parameter :: month = 1_IK
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
        real(RKG) :: julianDay
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
        integer(IKG) :: values(8)
        call date_and_time(values = values)
        CHECK_ASSERTION(__LINE__, all(values /= -huge(0_IKG)), SK_"@getJulianDay(): The processor does not have a clock.") ! fpp
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
        integer(IKG), parameter :: millisecond = 0_IK
#if     !getJulianDayYMDZHMS_ENABLED
        integer(IKG), parameter :: second = 0_IK
#if     !getJulianDayYMDZHM_ENABLED
        integer(IKG), parameter :: minute = 0_IK
#if     !getJulianDayYMDZH_ENABLED
        integer(IKG), parameter :: hour = 0_IK
#if     !getJulianDayYMDZ_ENABLED
        integer(IKG), parameter :: zone = 0_IK
#if     !getJulianDayYMD_ENABLED
        integer(IKG), parameter :: day = 1_IK
#if     !getJulianDayYM_ENABLED
        integer(IKG), parameter :: month = 1_IK
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
        real(RKG)   , parameter :: DAY_PER_HOUR = 1._RKG / 24._RKG
        real(RKG)   , parameter :: DAY_PER_MIN = 1._RKG / 1440._RKG
        real(RKG)   , parameter :: DAY_PER_SEC = 1._RKG / 86400._RKG
        integer(IKG), parameter :: VECTOR(12) = [306_IKG, 337_IKG, 0_IKG, 31_IKG, 61_IKG, 92_IKG, 122_IKG, 153_IKG, 184_IKG, 214_IKG, 245_IKG, 275_IKG]
        integer(IKG) :: yearCorrected, DateTimeUTC(8)

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

        if (DateTimeUTC(2) == 1_IKG .or. DateTimeUTC(2) == 2_IKG) then
            yearCorrected = DateTimeUTC(1) - 1_IKG
        else
            yearCorrected = DateTimeUTC(1)
        end if
        julianDay   = 1721118.5_RKG & ! LCOV_EXCL_LINE
                    + real(DateTimeUTC(3), RKG) & ! LCOV_EXCL_LINE
                    + real(VECTOR(DateTimeUTC(2)), RKG) & ! LCOV_EXCL_LINE
                    + 365._RKG * real(yearCorrected, RKG) & ! LCOV_EXCL_LINE
                    + real(floor(yearCorrected * 0.25_RKG, IKG) - floor(yearCorrected * 0.01_RKG, IKG) + floor(yearCorrected * 0.0025_RKG, IKG), RKG) & ! LCOV_EXCL_LINE
                    + DateTimeUTC(5) * DAY_PER_HOUR + DateTimeUTC(6) * DAY_PER_MIN + (real(DateTimeUTC(7),RKG) + 0.001_RKG * DateTimeUTC(8)) * DAY_PER_SEC
#endif

        !%%%%%%%%%%%%%%%%%%%%%%
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%
