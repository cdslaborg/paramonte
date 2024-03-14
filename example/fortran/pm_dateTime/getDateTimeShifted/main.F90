program example

    use pm_kind, only: IK, RK, SK
    use pm_io, only: display_type
    use pm_dateTime, only: getDateTimeShifted
    use pm_dateTime, only: getMillisecond
    use pm_dateTime, only: getDateTime
    use pm_dateTime, only: getSecond
    use pm_dateTime, only: getMinute
    use pm_dateTime, only: getMonth
    use pm_dateTime, only: getYear
    use pm_dateTime, only: getHour
    use pm_dateTime, only: getZone
    use pm_dateTime, only: getDay

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getDateTime()")
    call disp%show( getDateTime() )
    call disp%show("getDateTimeShifted(1._RK)")
    call disp%show( getDateTimeShifted(1._RK) )
    call disp%show("getDateTimeShifted(1._RK, getDateTime())")
    call disp%show( getDateTimeShifted(1._RK, getDateTime()) )

    call disp%skip()
    call disp%show("getDateTimeShifted(2._RK, 2000_IK)")
    call disp%show( getDateTimeShifted(2._RK, 2000_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(3._RK, 2000_IK, 1_IK)")
    call disp%show( getDateTimeShifted(3._RK, 2000_IK, 1_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(300._RK, 2000_IK, 1_IK, 1_IK)")
    call disp%show( getDateTimeShifted(300._RK, 2000_IK, 1_IK, 1_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(366._RK, 2000_IK, 1_IK, 1_IK, getZone())")
    call disp%show( getDateTimeShifted(366._RK, 2000_IK, 1_IK, 1_IK, getZone()) )

    call disp%skip()
    call disp%show("getDateTimeShifted(-366._RK, 2000_IK, 1_IK, 1_IK, getZone())")
    call disp%show( getDateTimeShifted(-366._RK, 2000_IK, 1_IK, 1_IK, getZone()) )

    call disp%skip()
    call disp%show("getDateTimeShifted(-366.25_RK, 2000_IK, 1_IK, 1_IK, 0_IK)")
    call disp%show( getDateTimeShifted(-366.25_RK, 2000_IK, 1_IK, 1_IK, 0_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(1._RK, getYear(), getMonth(), getDay(), getZone(), getHour())")
    call disp%show( getDateTimeShifted(1._RK, getYear(), getMonth(), getDay(), getZone(), getHour()) )

    call disp%skip()
    call disp%show("getDateTimeShifted(1._RK, getYear(), getMonth(), getDay(), getZone(), getHour(), getMinute())")
    call disp%show( getDateTimeShifted(1._RK, getYear(), getMonth(), getDay(), getZone(), getHour(), getMinute()) )

    call disp%skip()
    call disp%show("getDateTimeShifted(1._RK, getYear(), getMonth(), getDay(), getZone(), getHour(), getMinute(), getSecond())")
    call disp%show( getDateTimeShifted(1._RK, getYear(), getMonth(), getDay(), getZone(), getHour(), getMinute(), getSecond()) )

    call disp%skip()
    call disp%show("getDateTimeShifted(1._RK, getYear(), getMonth(), getDay(), getZone(), getHour(), getMinute(), getSecond(), getMillisecond())")
    call disp%show( getDateTimeShifted(1._RK, getYear(), getMonth(), getDay(), getZone(), getHour(), getMinute(), getSecond(), getMillisecond()) )

    call disp%skip()
    call disp%show("getDateTimeShifted(-1._RK, 2000_IK, 12_IK, 31_IK, -660_IK)")
    call disp%show( getDateTimeShifted(-1._RK, 2000_IK, 12_IK, 31_IK, -660_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(0._RK, 2000_IK, 12_IK, 31_IK, -660_IK)")
    call disp%show( getDateTimeShifted(0._RK, 2000_IK, 12_IK, 31_IK, -660_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(2*365._RK, 2000_IK, 12_IK, 31_IK, -660_IK)")
    call disp%show( getDateTimeShifted(2*365._RK, 2000_IK, 12_IK, 31_IK, -660_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTimeShifted(-2*365._RK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK)")
    call disp%show( getDateTimeShifted(-2*365._RK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(0._RK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK)")
    call disp%show( getDateTimeShifted(0._RK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(+2*365._RK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK)")
    call disp%show( getDateTimeShifted(+2*365._RK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTimeShifted(300._RK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK)")
    call disp%show( getDateTimeShifted(300._RK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(313._RK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK)")
    call disp%show( getDateTimeShifted(313._RK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(300._RK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK, 35_IK)")
    call disp%show( getDateTimeShifted(300._RK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK, 35_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(300._RK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK, 35_IK, 847_IK)")
    call disp%show( getDateTimeShifted(300._RK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK, 35_IK, 847_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(300._RK, 2000_IK, 2_IK, 29_IK, -660_IK, 18_IK, 21_IK, 35_IK, 847_IK)")
    call disp%show( getDateTimeShifted(300._RK, 2000_IK, 2_IK, 29_IK, -660_IK, 18_IK, 21_IK, 35_IK, 847_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(-300._RK, 2000_IK, 3_IK, 1_IK, +660_IK, 8_IK, 21_IK, 35_IK, 847_IK)")
    call disp%show( getDateTimeShifted(-300._RK, 2000_IK, 3_IK, 1_IK, +660_IK, 8_IK, 21_IK, 35_IK, 847_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(-100._RK, 1999_IK, 3_IK, 1_IK, +660_IK, 8_IK, 21_IK, 35_IK, 847_IK)")
    call disp%show( getDateTimeShifted(-100._RK, 1999_IK, 3_IK, 1_IK, +660_IK, 8_IK, 21_IK, 35_IK, 847_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(365._RK, 1999_IK, 1_IK, 1_IK, +660_IK, 0_IK, 0_IK, 0_IK, 0_IK)")
    call disp%show( getDateTimeShifted(365._RK, 1999_IK, 1_IK, 1_IK, +660_IK, 0_IK, 0_IK, 0_IK, 0_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(365._RK, 1999_IK, 1_IK, 1_IK, +660_IK)")
    call disp%show( getDateTimeShifted(365._RK, 1999_IK, 1_IK, 1_IK, +660_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(365._RK, 1999_IK, 1_IK, 1_IK)")
    call disp%show( getDateTimeShifted(365._RK, 1999_IK, 1_IK, 1_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(365._RK, 1999_IK, 1_IK)")
    call disp%show( getDateTimeShifted(365._RK, 1999_IK, 1_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(365._RK, 1999_IK)")
    call disp%show( getDateTimeShifted(365._RK, 1999_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(365._RK, 2000_IK)")
    call disp%show( getDateTimeShifted(365._RK, 2000_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(366._RK, 2000_IK)")
    call disp%show( getDateTimeShifted(366._RK, 2000_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(-366._RK, 2000_IK, 12_IK, 31_IK)")
    call disp%show( getDateTimeShifted(-366._RK, 2000_IK, 12_IK, 31_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(2.5_RK, -1_IK, 12_IK, 31_IK, -660_IK, 20_IK)")
    call disp%show( getDateTimeShifted(2.5_RK, -1_IK, 12_IK, 31_IK, -660_IK, 20_IK) )

    call disp%skip()
    call disp%show("getDateTimeShifted(-1.2_RK, 1_IK, 1_IK, 1_IK, +660_IK)")
    call disp%show( getDateTimeShifted(-1.2_RK, 1_IK, 1_IK, 1_IK, +660_IK) )

end program example