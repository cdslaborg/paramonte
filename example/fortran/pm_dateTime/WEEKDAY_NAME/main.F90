program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: WEEKDAY_NAME, getWeekDay

    implicit none

    integer :: i
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("WEEKDAY_NAME")
    call disp%show( WEEKDAY_NAME , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekDay()")
    call disp%show( getWeekDay() )
    call disp%show("trim(WEEKDAY_NAME(getWeekDay())) ! current day of week.")
    call disp%show( trim(WEEKDAY_NAME(getWeekDay())) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekDay()")
    call disp%show( getWeekDay() )
    call disp%show("WEEKDAY_NAME(getWeekDay())(1:3)")
    call disp%show( WEEKDAY_NAME(getWeekDay())(1:3) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("[( WEEKDAY_NAME(i)(1:3), i = 0, 6 )]")
    call disp%show( [( WEEKDAY_NAME(i)(1:3), i = 0, 6 )] , deliml = SK_"""" )
    call disp%skip()

end program example