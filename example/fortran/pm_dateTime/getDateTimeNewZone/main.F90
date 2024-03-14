program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: getDateTime
    use pm_dateTime, only: getDateTimeUTC
    use pm_dateTime, only: getDateTimeNewZone

    implicit none

    integer(IK) :: DateTime(8)
    integer(IK) :: DateTimeUTC(8)
    integer(IK) :: DateTimeNewZone(8)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("DateTime = getDateTime()")
                    DateTime = getDateTime()
    call disp%show("DateTime")
    call disp%show( DateTime )
    call disp%show("DateTimeNewZone = getDateTimeNewZone(newzone = DateTime(4))")
                    DateTimeNewZone = getDateTimeNewZone(newzone = DateTime(4))
    call disp%show("DateTimeNewZone")
    call disp%show( DateTimeNewZone )
    call disp%skip()

    call disp%skip()
    call disp%show("DateTime = getDateTime()")
                    DateTime = getDateTime()
    call disp%show("DateTime")
    call disp%show( DateTime )
    call disp%show("DateTimeNewZone = getDateTimeNewZone(newzone = -DateTime(4))")
                    DateTimeNewZone = getDateTimeNewZone(newzone = -DateTime(4))
    call disp%show("DateTimeNewZone")
    call disp%show( DateTimeNewZone )
    call disp%skip()

    call disp%skip()
    call disp%show("DateTimeUTC = getDateTimeUTC()")
                    DateTimeUTC = getDateTimeUTC()
    call disp%show("DateTimeUTC")
    call disp%show( DateTimeUTC )
    call disp%show("DateTimeNewZone = getDateTimeNewZone(newzone = DateTimeUTC(4), Values = DateTimeUTC)")
                    DateTimeNewZone = getDateTimeNewZone(newzone = DateTimeUTC(4), Values = DateTimeUTC)
    call disp%show("DateTimeNewZone")
    call disp%show( DateTimeNewZone )
    call disp%skip()

    call disp%skip()
    call disp%show("DateTimeUTC = getDateTimeUTC()")
                    DateTimeUTC = getDateTimeUTC()
    call disp%show("DateTimeUTC")
    call disp%show( DateTimeUTC )
    call disp%show("DateTimeNewZone = getDateTimeNewZone(newzone = 0_IK, Values = DateTimeUTC)")
                    DateTimeNewZone = getDateTimeNewZone(newzone = 0_IK, Values = DateTimeUTC)
    call disp%show("DateTimeNewZone")
    call disp%show( DateTimeNewZone )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTimeNewZone(-660_IK, 2000_IK, 12_IK, 31_IK, -660_IK)")
    call disp%show( getDateTimeNewZone(-660_IK, 2000_IK, 12_IK, 31_IK, -660_IK) )

    call disp%skip()
    call disp%show("getDateTimeNewZone(0_IK, 2000_IK, 12_IK, 31_IK, -660_IK)")
    call disp%show( getDateTimeNewZone(0_IK, 2000_IK, 12_IK, 31_IK, -660_IK) )

    call disp%skip()
    call disp%show("getDateTimeNewZone(+660_IK, 2000_IK, 12_IK, 31_IK, -660_IK)")
    call disp%show( getDateTimeNewZone(+660_IK, 2000_IK, 12_IK, 31_IK, -660_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTimeNewZone(-660_IK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK)")
    call disp%show( getDateTimeNewZone(-660_IK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK) )

    call disp%skip()
    call disp%show("getDateTimeNewZone(0_IK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK)")
    call disp%show( getDateTimeNewZone(0_IK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK) )

    call disp%skip()
    call disp%show("getDateTimeNewZone(+660_IK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK)")
    call disp%show( getDateTimeNewZone(+660_IK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTimeNewZone(300_IK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK)")
    call disp%show( getDateTimeNewZone(300_IK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK) )

    call disp%skip()
    call disp%show("getDateTimeNewZone(313_IK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK)")
    call disp%show( getDateTimeNewZone(313_IK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK) )

    call disp%skip()
    call disp%show("getDateTimeNewZone(300_IK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK, 35_IK)")
    call disp%show( getDateTimeNewZone(300_IK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK, 35_IK) )

    call disp%skip()
    call disp%show("getDateTimeNewZone(300_IK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK, 35_IK, 847_IK)")
    call disp%show( getDateTimeNewZone(300_IK, 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK, 35_IK, 847_IK) )

    call disp%skip()
    call disp%show("getDateTimeNewZone(300_IK, 2000_IK, 2_IK, 29_IK, -660_IK, 18_IK, 21_IK, 35_IK, 847_IK)")
    call disp%show( getDateTimeNewZone(300_IK, 2000_IK, 2_IK, 29_IK, -660_IK, 18_IK, 21_IK, 35_IK, 847_IK) )

    call disp%skip()
    call disp%show("getDateTimeNewZone(-300_IK, 2000_IK, 3_IK, 1_IK, +660_IK, 8_IK, 21_IK, 35_IK, 847_IK)")
    call disp%show( getDateTimeNewZone(-300_IK, 2000_IK, 3_IK, 1_IK, +660_IK, 8_IK, 21_IK, 35_IK, 847_IK) )

    call disp%skip()
    call disp%show("getDateTimeNewZone(-300_IK, 1999_IK, 3_IK, 1_IK, +660_IK, 8_IK, 21_IK, 35_IK, 847_IK)")
    call disp%show( getDateTimeNewZone(-300_IK, 1999_IK, 3_IK, 1_IK, +660_IK, 8_IK, 21_IK, 35_IK, 847_IK) )

    call disp%skip()
    call disp%show("getDateTimeNewZone(300_IK, -1_IK, 12_IK, 31_IK, -660_IK, 20_IK)")
    call disp%show( getDateTimeNewZone(300_IK, -1_IK, 12_IK, 31_IK, -660_IK, 20_IK) )

    call disp%skip()
    call disp%show("getDateTimeNewZone(300_IK, 1_IK, 1_IK, 1_IK, +660_IK)")
    call disp%show( getDateTimeNewZone(300_IK, 1_IK, 1_IK, 1_IK, +660_IK) )

end program example