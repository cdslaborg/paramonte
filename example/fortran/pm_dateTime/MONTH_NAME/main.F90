program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: MONTH_NAME, getMonth

    implicit none

    integer :: i
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("MONTH_NAME")
    call disp%show( MONTH_NAME , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getMonth()")
    call disp%show( getMonth() )
    call disp%show("trim(MONTH_NAME(getMonth()))")
    call disp%show( trim(MONTH_NAME(getMonth())) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getMonth()")
    call disp%show( getMonth() )
    call disp%show("MONTH_NAME(getMonth())(1:3)")
    call disp%show( MONTH_NAME(getMonth())(1:3) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("[( MONTH_NAME(i)(1:3), i = 1, 12 )]")
    call disp%show( [( MONTH_NAME(i)(1:3), i = 1, 12 )] , deliml = SK_"""" )
    call disp%skip()

end program example