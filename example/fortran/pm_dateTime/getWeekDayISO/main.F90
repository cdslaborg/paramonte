program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: WEEKDAY_NAME_ISO
    use pm_dateTime, only: getWeekDayISO

    implicit none

    integer(IK) :: Values(8)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Get the day number of the week of a Gregorian calendar date, assuming Monday is the first day of the week.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekDayISO()")
    call disp%show( getWeekDayISO() )
    call disp%show("trim(WEEKDAY_NAME_ISO(getWeekDayISO())) ! current day of week.")
    call disp%show( trim(WEEKDAY_NAME_ISO(getWeekDayISO())) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call date_and_time(values = Values(1:8))")
                    call date_and_time(values = Values(1:8))
    call disp%show("getWeekDayISO(Values(1:3))")
    call disp%show( getWeekDayISO(Values(1:3)) )
    call disp%show("getWeekDayISO(year = Values(1), month = Values(2), day = Values(3))")
    call disp%show( getWeekDayISO(year = Values(1), month = Values(2), day = Values(3)) )
    call disp%show("SK_'Today is '//trim(WEEKDAY_NAME_ISO(getWeekDayISO(Values(1:3))))//SK_'.'")
    call disp%show( SK_'Today is '//trim(WEEKDAY_NAME_ISO(getWeekDayISO(Values(1:3))))//SK_'.' , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekDayISO(1582_IK, 10_IK, 15_IK)")
    call disp%show( getWeekDayISO(1582_IK, 10_IK, 15_IK) )
    call disp%show("trim(WEEKDAY_NAME_ISO(getWeekDayISO(1582_IK, 10_IK, 15_IK)))//SK_' October 15, 1582 is the first day of the Gregorian Calendar.'")
    call disp%show( trim(WEEKDAY_NAME_ISO(getWeekDayISO(1582_IK, 10_IK, 15_IK)))//SK_' October 15, 1582 is the first day of the Gregorian Calendar.' , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekDayISO(1_IK, 1_IK, 1_IK)")
    call disp%show( getWeekDayISO(1_IK, 1_IK, 1_IK) )
    call disp%show("SK_'January 1, 0001 is a '//trim(WEEKDAY_NAME_ISO(getWeekDayISO(1_IK, 1_IK, 1_IK)))//SK_' in the proleptic Gregorian Calendar.'")
    call disp%show( SK_'January 1, 0001 is a '//trim(WEEKDAY_NAME_ISO(getWeekDayISO(1_IK, 1_IK, 1_IK)))//SK_' in the proleptic Gregorian Calendar.' , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekDayISO(2000_IK, 1_IK, 1_IK)")
    call disp%show( getWeekDayISO(2000_IK, 1_IK, 1_IK) )
    call disp%show("trim(WEEKDAY_NAME_ISO(getWeekDayISO(2000_IK, 1_IK, 1_IK)))")
    call disp%show( trim(WEEKDAY_NAME_ISO(getWeekDayISO(2000_IK, 1_IK, 1_IK))) , deliml = SK_"""" )
    call disp%skip()

end program example