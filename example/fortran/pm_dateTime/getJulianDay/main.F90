program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: getJulianDay

    implicit none

    integer(IK) :: Values(8)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getJulianDay()")
    call disp%show( getJulianDay() )
    call disp%skip()

    call disp%skip()
    call disp%show("call date_and_time(values = Values)")
                    call date_and_time(values = Values)
    call disp%show("getJulianDay(Values(1))")
    call disp%show( getJulianDay(Values(1)) )
    call disp%show("getJulianDay(Values(1:2))")
    call disp%show( getJulianDay(Values(1:2)) )
    call disp%show("getJulianDay(Values(1:3))")
    call disp%show( getJulianDay(Values(1:3)) )
    call disp%show("getJulianDay(Values(1:4))")
    call disp%show( getJulianDay(Values(1:4)) )
    call disp%show("getJulianDay(Values(1:5))")
    call disp%show( getJulianDay(Values(1:5)) )
    call disp%show("getJulianDay(Values(1:6))")
    call disp%show( getJulianDay(Values(1:6)) )
    call disp%show("getJulianDay(Values(1:7))")
    call disp%show( getJulianDay(Values(1:7)) )
    call disp%show("getJulianDay(Values(1:8))")
    call disp%show( getJulianDay(Values(1:8)) )
    call disp%skip()

    call disp%skip()
    call disp%show("getJulianDay(-4713_IK, 11_IK, 24_IK) ! -0.5")
    call disp%show( getJulianDay(-4713_IK, 11_IK, 24_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getJulianDay(-4713_IK, 11_IK, 24_IK, 0_IK, 12_IK) ! 0.0 (until noon)")
    call disp%show( getJulianDay(-4713_IK, 11_IK, 24_IK, 0_IK, 12_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getJulianDay(-4713_IK, 11_IK, 25_IK, 0_IK, 0_IK) ! 0.5")
    call disp%show( getJulianDay(-4713_IK, 11_IK, 25_IK, 0_IK, 0_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getJulianDay(-1_IK, 1_IK, 1_IK) ! 1720694.5")
    call disp%show( getJulianDay(-1_IK, 1_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getJulianDay(1_IK, 1_IK, 1_IK) ! 1721425.5")
    call disp%show( getJulianDay(1_IK, 1_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getJulianDay(1582_IK, 10_IK, 15_IK) ! 2299160.5 (The first day of the Gregorian Calendar reform)")
    call disp%show( getJulianDay(1582_IK, 10_IK, 15_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getJulianDay(1901_IK, 1_IK, 1_IK) ! 2415385.5 (The start of the 20th century)")
    call disp%show( getJulianDay(1901_IK, 1_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getJulianDay(1970_IK, 1_IK, 1_IK) ! 2440587.5 (The Unix reference date)")
    call disp%show( getJulianDay(1970_IK, 1_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getJulianDay(2000_IK, 2_IK, 28_IK) ! 2451602.5")
    call disp%show( getJulianDay(2000_IK, 2_IK, 28_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getJulianDay(2000_IK, 1_IK, 1_IK) ! 2451544.5")
    call disp%show( getJulianDay(2000_IK, 1_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getJulianDay(2000_IK, 1_IK, 1_IK, 0_IK, 18_IK) ! 2451545.25")
    call disp%show( getJulianDay(2000_IK, 1_IK, 1_IK, 0_IK, 18_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getJulianDay(2022_IK, 5_IK, 8_IK) ! 2459707.5")
    call disp%show( getJulianDay(2022_IK, 5_IK, 8_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getJulianDay(3000_IK, 12_IK, 31_IK) ! 2817151.5")
    call disp%show( getJulianDay(3000_IK, 12_IK, 31_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getJulianDay(9999_IK, 12_IK, 31_IK) ! 5373483.5")
    call disp%show( getJulianDay(9999_IK, 12_IK, 31_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getJulianDay(99999_IK, 12_IK, 31_IK) ! 38245308.5")
    call disp%show( getJulianDay(99999_IK, 12_IK, 31_IK) )
    call disp%skip()

end program example