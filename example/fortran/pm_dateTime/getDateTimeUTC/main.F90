program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: getDateTime
    use pm_dateTime, only: getDateTimeUTC

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getDateTime()")
    call disp%show( getDateTime() )
    call disp%show("getDateTimeUTC()")
    call disp%show( getDateTimeUTC() )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTimeUTC(2000_IK, 12_IK, 31_IK, -660_IK)")
    call disp%show( getDateTimeUTC(2000_IK, 12_IK, 31_IK, -660_IK) )

    call disp%skip()
    call disp%show("getDateTimeUTC(2000_IK, 12_IK, 31_IK, -660_IK, 18_IK)")
    call disp%show( getDateTimeUTC(2000_IK, 12_IK, 31_IK, -660_IK, 18_IK) )

    call disp%skip()
    call disp%show("getDateTimeUTC(2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK)")
    call disp%show( getDateTimeUTC(2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK) )

    call disp%skip()
    call disp%show("getDateTimeUTC(2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK, 35_IK)")
    call disp%show( getDateTimeUTC(2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK, 35_IK) )

    call disp%skip()
    call disp%show("getDateTimeUTC(2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK, 35_IK, 847_IK)")
    call disp%show( getDateTimeUTC(2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK, 35_IK, 847_IK) )

    call disp%skip()
    call disp%show("getDateTimeUTC(2000_IK, 2_IK, 29_IK, -660_IK, 18_IK, 21_IK, 35_IK, 847_IK)")
    call disp%show( getDateTimeUTC(2000_IK, 2_IK, 29_IK, -660_IK, 18_IK, 21_IK, 35_IK, 847_IK) )

    call disp%skip()
    call disp%show("getDateTimeUTC(2000_IK, 3_IK, 1_IK, +660_IK, 8_IK, 21_IK, 35_IK, 847_IK)")
    call disp%show( getDateTimeUTC(2000_IK, 3_IK, 1_IK, +660_IK, 8_IK, 21_IK, 35_IK, 847_IK) )

    call disp%skip()
    call disp%show("getDateTimeUTC(1999_IK, 3_IK, 1_IK, +660_IK, 8_IK, 21_IK, 35_IK, 847_IK)")
    call disp%show( getDateTimeUTC(1999_IK, 3_IK, 1_IK, +660_IK, 8_IK, 21_IK, 35_IK, 847_IK) )

    call disp%skip()
    call disp%show("getDateTimeUTC(-1_IK, 12_IK, 31_IK, -660_IK, 20_IK)")
    call disp%show( getDateTimeUTC(-1_IK, 12_IK, 31_IK, -660_IK, 20_IK) )

    call disp%skip()
    call disp%show("getDateTimeUTC(1_IK, 1_IK, 1_IK, +660_IK)")
    call disp%show( getDateTimeUTC(1_IK, 1_IK, 1_IK, +660_IK) )

end program example