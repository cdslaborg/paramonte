program example

    use pm_kind, only: IK, RK, SK
    use pm_io, only: display_type
    use pm_dateTime, only: getDateTime
    use pm_dateTime, only: isMorning

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getDateTime()")
    call disp%show( getDateTime() )
    call disp%show("isMorning() ! is it morning at local time?")
    call disp%show( isMorning() )
    call disp%skip()

    call disp%skip()
    call disp%show("isMorning(0) ! is it morning at UTC?")
    call disp%show( isMorning(0) )
    call disp%skip()

    call disp%skip()
    call disp%show("isMorning(60 * [-12, 3, 14]) ! is it morning in the specified timezones with respect to UTC?")
    call disp%show( isMorning(60 * [-12, 3, 14]) )
    call disp%skip()

    call disp%skip()
    call disp%show("isMorning(-0.5_RK) ! November 24, -4713, start of the day (just after midnight)")
    call disp%show( isMorning(-0.5_RK) )
    call disp%skip()

    call disp%skip()
    call disp%show("isMorning(-0.3_RK) ! November 24, -4713, start of the day (after midnight, before noon)")
    call disp%show( isMorning(-0.3_RK) )
    call disp%skip()

    call disp%skip()
    call disp%show("isMorning(0._RK) ! November 24, -4713, noon (start of JD in Gregorian Calendar)")
    call disp%show( isMorning(0._RK) )
    call disp%skip()

    call disp%skip()
    call disp%show("isMorning(0.5_RK) ! November 25, -4713, start of the day (just after midnight)")
    call disp%show( isMorning(0.5_RK) )
    call disp%skip()

    call disp%skip()
    call disp%show("isMorning(0.5_RK, zone = -6 * 60) ! November 25, -4713, start of the UTC day (just after midnight) but at timezone -6 (Texas) where it is still afternoon")
    call disp%show( isMorning(0.5_RK, zone = -6 * 60) )
    call disp%skip()

    call disp%skip()
    call disp%show("isMorning([2440587.5_RK, 2440587.9_RK, 2444239.49999_RK]) ! morning immediately after midnight, morning before noon, afternoon right before midnight")
    call disp%show( isMorning([2440587.5_RK, 2440587.9_RK, 2444239.49999_RK]) )
    call disp%skip()

    call disp%skip()
    call disp%show("isMorning([2440587.5_RK, 2440587.9_RK, 2444239.49999_RK], zone = 60 * [-1, +14, -12]) ! afternoon near midnight, afternoon the day before UTC, morning right before none at time zone -720")
    call disp%show( isMorning([2440587.5_RK, 2440587.9_RK, 2444239.49999_RK], zone = 60 * [-1, +14, -12]) )
    call disp%skip()

end program example