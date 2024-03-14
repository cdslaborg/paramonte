program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: WEEKDAY_NAME_ISO, getWeekDayISO

    implicit none

    integer :: i
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("WEEKDAY_NAME_ISO")
    call disp%show( WEEKDAY_NAME_ISO , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekDayISO()")
    call disp%show( getWeekDayISO() )
    call disp%show("trim(WEEKDAY_NAME_ISO(getWeekDayISO())) ! current day of week.")
    call disp%show( trim(WEEKDAY_NAME_ISO(getWeekDayISO())) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekDayISO()")
    call disp%show( getWeekDayISO() )
    call disp%show("WEEKDAY_NAME_ISO(getWeekDayISO())(1:3)")
    call disp%show( WEEKDAY_NAME_ISO(getWeekDayISO())(1:3) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("[( WEEKDAY_NAME_ISO(i)(1:3), i = 1, 7 )]")
    call disp%show( [( WEEKDAY_NAME_ISO(i)(1:3), i = 1, 7 )] , deliml = SK_"""" )
    call disp%skip()

end program example