program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: getDateTime
    use pm_dateTime, only: getOrdinalDay

    implicit none

    integer(IK) :: Values(8)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getDateTime()")
    call disp%show( getDateTime() )
    call disp%show("getOrdinalDay()")
    call disp%show( getOrdinalDay() )
    call disp%skip()

    call disp%skip()
    call disp%show("call date_and_time(values = Values)")
                    call date_and_time(values = Values)
    call disp%show("Values(1:3)")
    call disp%show( Values(1:3) )
    call disp%show("getOrdinalDay(Values(1:3))")
    call disp%show( getOrdinalDay(Values(1:3)) )
    call disp%skip()

    call disp%skip()
    call disp%show("getOrdinalDay(-4713_IK, 11_IK, 24_IK)")
    call disp%show( getOrdinalDay(-4713_IK, 11_IK, 24_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getOrdinalDay(-1_IK, 1_IK, 1_IK)")
    call disp%show( getOrdinalDay(-1_IK, 1_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getOrdinalDay(-1_IK, 12_IK, 31_IK)")
    call disp%show( getOrdinalDay(-1_IK, 12_IK, 31_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getOrdinalDay(0_IK, 12_IK, 31_IK)")
    call disp%show( getOrdinalDay(0_IK, 12_IK, 31_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getOrdinalDay(1_IK, 12_IK, 31_IK)")
    call disp%show( getOrdinalDay(1_IK, 12_IK, 31_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getOrdinalDay(1582_IK, 10_IK, 15_IK)")
    call disp%show( getOrdinalDay(1582_IK, 10_IK, 15_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getOrdinalDay(1901_IK, 1_IK, 1_IK)")
    call disp%show( getOrdinalDay(1901_IK, 1_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getOrdinalDay(1999_IK, 3_IK, 1_IK)")
    call disp%show( getOrdinalDay(1999_IK, 3_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getOrdinalDay(2000_IK, 3_IK, 1_IK)")
    call disp%show( getOrdinalDay(2000_IK, 3_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getOrdinalDay(1999_IK, 4_IK, 15_IK) ! 105")
    call disp%show( getOrdinalDay(1999_IK, 4_IK, 15_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getOrdinalDay(2000_IK, 4_IK, 15_IK) ! 106 (leap year)")
    call disp%show( getOrdinalDay(2000_IK, 4_IK, 15_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getOrdinalDay(9999_IK, 12_IK, 31_IK)")
    call disp%show( getOrdinalDay(9999_IK, 12_IK, 31_IK) )
    call disp%skip()

end program example