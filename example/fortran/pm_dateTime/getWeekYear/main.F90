program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: getDateTime
    use pm_dateTime, only: getWeekYear

    implicit none

    integer(IK) :: Values(8)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getDateTime()")
    call disp%show( getDateTime() )
    call disp%show("getWeekYear()")
    call disp%show( getWeekYear() )
    call disp%skip()

    call disp%skip()
    call disp%show("call date_and_time(values = Values)")
                    call date_and_time(values = Values)
    call disp%show("Values(1:3)")
    call disp%show( Values(1:3) )
    call disp%show("getWeekYear(Values(1:3))")
    call disp%show( getWeekYear(Values(1:3)) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekYear(1977_IK, 1_IK, 1_IK) ! 1976-W53-6")
    call disp%show( getWeekYear(1977_IK, 1_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekYear(1977_IK, 1_IK, 2_IK) ! 1976-W53-7")
    call disp%show( getWeekYear(1977_IK, 1_IK, 2_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekYear(1977_IK, 12_IK, 31_IK) ! 1977-W52-6")
    call disp%show( getWeekYear(1977_IK, 12_IK, 31_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekYear(1978_IK, 1_IK, 1_IK) ! 1977-W52-7")
    call disp%show( getWeekYear(1978_IK, 1_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekYear(1978_IK, 1_IK, 2_IK) ! 1978-W01-1")
    call disp%show( getWeekYear(1978_IK, 1_IK, 2_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekYear(1978_IK, 12_IK, 31_IK) ! 1978-W52-7")
    call disp%show( getWeekYear(1978_IK, 12_IK, 31_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekYear(1979_IK, 1_IK, 1_IK) ! 1979-W01-1")
    call disp%show( getWeekYear(1979_IK, 1_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekYear(1979_IK, 12_IK, 30_IK) ! 1979-W52-7")
    call disp%show( getWeekYear(1979_IK, 12_IK, 30_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getWeekYear(1979_IK, 12_IK, 31_IK) ! 1980-W01-1")
    call disp%show( getWeekYear(1979_IK, 12_IK, 31_IK) )
    call disp%skip()

end program example