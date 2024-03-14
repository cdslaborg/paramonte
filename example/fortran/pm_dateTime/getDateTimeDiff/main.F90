program example

    use pm_kind, only: IK, RK, SK
    use pm_io, only: display_type
    use pm_dateTime, only: getDateTime
    use pm_dateTime, only: getDateTimeUTC
    use pm_dateTime, only: getDateTimeDiff
    use pm_dateTime, only: getDateTimeShifted

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getDateTimeUTC()")
    call disp%show( getDateTimeUTC() )
    call disp%show("getDateTime(2000_IK)")
    call disp%show( getDateTime(2000_IK) )
    call disp%show("getDateTimeDiff(getDateTime(), getDateTime(2000_IK))")
    call disp%show( getDateTimeDiff(getDateTime(), getDateTime(2000_IK)) )
    call disp%show("getDateTimeDiff(getDateTime(), getDateTimeShifted(-1.5_RK))")
    call disp%show( getDateTimeDiff(getDateTime(), getDateTimeShifted(-1.5_RK)) )

    call disp%skip()
    call disp%show("getDateTime(1999_IK, 2_IK, 28_IK)")
    call disp%show( getDateTime(1999_IK, 2_IK, 28_IK) )
    call disp%show("getDateTime(1999_IK, 3_IK, 1_IK)")
    call disp%show( getDateTime(1999_IK, 3_IK, 1_IK) )
    call disp%show("getDateTimeDiff(getDateTime(1999_IK, 3_IK, 1_IK), getDateTime(1999_IK, 2_IK, 28_IK)) ! year 1999 is not a leap year.")
    call disp%show( getDateTimeDiff(getDateTime(1999_IK, 3_IK, 1_IK), getDateTime(1999_IK, 2_IK, 28_IK)) )

    call disp%skip()
    call disp%show("getDateTime(2000_IK, 2_IK, 28_IK)")
    call disp%show( getDateTime(2000_IK, 2_IK, 28_IK) )
    call disp%show("getDateTime(2000_IK, 3_IK, 1_IK)")
    call disp%show( getDateTime(2000_IK, 3_IK, 1_IK) )
    call disp%show("getDateTimeDiff(getDateTime(2000_IK, 3_IK, 1_IK), getDateTime(2000_IK, 2_IK, 28_IK)) ! year 2000 is a leap year.")
    call disp%show( getDateTimeDiff(getDateTime(2000_IK, 3_IK, 1_IK), getDateTime(2000_IK, 2_IK, 28_IK)) )

    call disp%skip()
    call disp%show("getDateTime(2019_IK, 1_IK, 1_IK)")
    call disp%show( getDateTime(2019_IK, 1_IK, 1_IK) )
    call disp%show("getDateTime(2019_IK, 12_IK, 31_IK)")
    call disp%show( getDateTime(2019_IK, 12_IK, 31_IK) )
    call disp%show("getDateTimeDiff(getDateTime(2019_IK, 1_IK, 1_IK), getDateTime(2019_IK, 12_IK, 31_IK)) ! year 2019 is not a leap year.")
    call disp%show( getDateTimeDiff(getDateTime(2019_IK, 1_IK, 1_IK), getDateTime(2019_IK, 12_IK, 31_IK)) )

    call disp%skip()
    call disp%show("getDateTime(2020_IK, 1_IK, 1_IK)")
    call disp%show( getDateTime(2020_IK, 1_IK, 1_IK) )
    call disp%show("getDateTime(2020_IK, 12_IK, 31_IK)")
    call disp%show( getDateTime(2020_IK, 12_IK, 31_IK) )
    call disp%show("getDateTimeDiff(getDateTime(2020_IK, 1_IK, 1_IK), getDateTime(2020_IK, 12_IK, 31_IK)) ! year 2020 is a leap year.")
    call disp%show( getDateTimeDiff(getDateTime(2020_IK, 1_IK, 1_IK), getDateTime(2020_IK, 12_IK, 31_IK)) )

    call disp%skip()
    call disp%show("getDateTime(2020_IK, 1_IK, 1_IK)")
    call disp%show( getDateTime(2020_IK, 1_IK, 1_IK) )
    call disp%show("getDateTime(2021_IK, 1_IK, 1_IK)")
    call disp%show( getDateTime(2021_IK, 1_IK, 1_IK) )
    call disp%show("getDateTimeDiff(getDateTime(2020_IK, 1_IK, 1_IK), getDateTime(2021_IK, 1_IK, 1_IK)) ! year 2020 is a leap year.")
    call disp%show( getDateTimeDiff(getDateTime(2020_IK, 1_IK, 1_IK), getDateTime(2021_IK, 1_IK, 1_IK)) )

    call disp%skip()
    call disp%show("getDateTime(1_IK, 1_IK, 1_IK, zone = +12_IK * 60_IK)")
    call disp%show( getDateTime(1_IK, 1_IK, 1_IK, zone = +12_IK * 60_IK) )
    call disp%show("getDateTime(1_IK, 1_IK, 1_IK, zone = -12_IK * 60_IK)")
    call disp%show( getDateTime(1_IK, 1_IK, 1_IK, zone = -12_IK * 60_IK) )
    call disp%show("getDateTimeDiff(getDateTime(1_IK, 1_IK, 1_IK, zone = +12_IK * 60_IK), getDateTime(1_IK, 1_IK, 1_IK, zone = -12_IK * 60_IK))")
    call disp%show( getDateTimeDiff(getDateTime(1_IK, 1_IK, 1_IK, zone = +12_IK * 60_IK), getDateTime(1_IK, 1_IK, 1_IK, zone = -12_IK * 60_IK)) )

    call disp%skip()
    call disp%show("getDateTime(-1_IK, 1_IK, 1_IK, zone = +12_IK * 60_IK)")
    call disp%show( getDateTime(-1_IK, 1_IK, 1_IK, zone = +12_IK * 60_IK) )
    call disp%show("getDateTime(-1_IK, 1_IK, 1_IK, zone = -12_IK * 60_IK)")
    call disp%show( getDateTime(-1_IK, 1_IK, 1_IK, zone = -12_IK * 60_IK) )
    call disp%show("getDateTimeDiff(getDateTime(-1_IK, 1_IK, 1_IK, zone = -12_IK * 60_IK), getDateTime(-1_IK, 1_IK, 1_IK, zone = +12_IK * 60_IK))")
    call disp%show( getDateTimeDiff(getDateTime(-1_IK, 1_IK, 1_IK, zone = -12_IK * 60_IK), getDateTime(-1_IK, 1_IK, 1_IK, zone = +12_IK * 60_IK)) )

end program example