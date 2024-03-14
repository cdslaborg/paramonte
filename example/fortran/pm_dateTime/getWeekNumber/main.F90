program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: getDateTime
    use pm_dateTime, only: getWeekNumber

    implicit none

    integer(IK) :: Values(8)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getDateTime()")
    call disp%show( getDateTime() )
    call disp%show("getWeekNumber()")
    call disp%show( getWeekNumber() )
    call disp%skip()

    call disp%skip()
    call disp%show("call date_and_time(values = Values)")
                    call date_and_time(values = Values)
    call disp%show("Values(1:3)")
    call disp%show( Values(1:3) )
    call disp%show("getWeekNumber(Values(1:3))")
    call disp%show( getWeekNumber(Values(1:3)) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekNumber(+2022, +5, +8) ! Sunday: 18 weeks so far with majority of days in 2022.")
    call disp%show( getWeekNumber(+2022, +5, +8) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekNumber(+2022, +5, +9) ! Monday: 19 weeks so far with majority of days in 2022.")
    call disp%show( getWeekNumber(+2022, +5, +9) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekNumber(+2022, +5, +12) ! Thursday: 19 weeks so far with majority of days in 2022.")
    call disp%show( getWeekNumber(+2022, +5, +12) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekNumber(-4713_IK, 11_IK, 24_IK)")
    call disp%show( getWeekNumber(-4713_IK, 11_IK, 24_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekNumber(-1_IK, 1_IK, 1_IK) ! The last week of 3 BC.")
    call disp%show( getWeekNumber(-1_IK, 1_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekNumber(-1_IK, 12_IK, 31_IK) ! The last week of 2 BC.")
    call disp%show( getWeekNumber(-1_IK, 12_IK, 31_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekNumber(0_IK, 12_IK, 31_IK) ! The last week of 1 BC.")
    call disp%show( getWeekNumber(0_IK, 12_IK, 31_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekNumber(1_IK, 12_IK, 31_IK) ! The first week of year 2 AD")
    call disp%show( getWeekNumber(1_IK, 12_IK, 31_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekNumber(1582_IK, 10_IK, 15_IK)")
    call disp%show( getWeekNumber(1582_IK, 10_IK, 15_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekNumber(1901_IK, 1_IK, 1_IK)")
    call disp%show( getWeekNumber(1901_IK, 1_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekNumber(1999_IK, 3_IK, 1_IK) ! 9 weeks until March 1 with majority of days in 1999.")
    call disp%show( getWeekNumber(1999_IK, 3_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekNumber(2000_IK, 3_IK, 1_IK)")
    call disp%show( getWeekNumber(2000_IK, 3_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekNumber(1999_IK, 4_IK, 15_IK)")
    call disp%show( getWeekNumber(1999_IK, 4_IK, 15_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekNumber(2000_IK, 4_IK, 15_IK)")
    call disp%show( getWeekNumber(2000_IK, 4_IK, 15_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekNumber(9999_IK, 12_IK, 31_IK)")
    call disp%show( getWeekNumber(9999_IK, 12_IK, 31_IK) )
    call disp%skip()

end program example