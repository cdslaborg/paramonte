program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: DAYS_OF_MONTH_LEAP

    implicit none

    integer :: i
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("DAYS_OF_MONTH_LEAP")
    call disp%show( DAYS_OF_MONTH_LEAP )
    call disp%skip()

    call disp%skip()
    call disp%show("sum(DAYS_OF_MONTH_LEAP)")
    call disp%show( sum(DAYS_OF_MONTH_LEAP) )
    call disp%skip()

end program example