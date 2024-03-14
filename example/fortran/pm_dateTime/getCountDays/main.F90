program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: getCountDays
    use pm_dateTime, only: isLeapYear
    use pm_dateTime, only: getMonth
    use pm_dateTime, only: getYear

    implicit none

    integer(IK), allocatable :: Year(:)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Get the count of days in a given year of the Gregorian calendar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getYear()")
    call disp%show( getYear() )
    call disp%show("isLeapYear(getYear())")
    call disp%show( isLeapYear(getYear()) )
    call disp%show("getCountDays(getYear())")
    call disp%show( getCountDays(getYear()) )
    call disp%skip()

    call disp%skip()
    call disp%show("isLeapYear(2000_IK)")
    call disp%show( isLeapYear(2000_IK) )
    call disp%show("getCountDays(2000_IK)")
    call disp%show( getCountDays(2000_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("isLeapYear(2024_IK)")
    call disp%show( isLeapYear(2024_IK) )
    call disp%show("getCountDays(2024_IK)")
    call disp%show( getCountDays(2024_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("Year = [1_IK, 2_IK, 3_IK, 4_IK]")
                    Year = [1_IK, 2_IK, 3_IK, 4_IK]
    call disp%show("Year")
    call disp%show( Year )
    call disp%show("isLeapYear(Year)")
    call disp%show( isLeapYear(Year) )
    call disp%show("getCountDays(Year)")
    call disp%show( getCountDays(Year) )
    call disp%skip()

    call disp%skip()
    call disp%show("Year = -[1_IK, 2_IK, 3_IK, 4_IK]")
                    Year = -[1_IK, 2_IK, 3_IK, 4_IK]
    call disp%show("Year")
    call disp%show( Year )
    call disp%show("isLeapYear(Year)")
    call disp%show( isLeapYear(Year) )
    call disp%show("getCountDays(Year)")
    call disp%show( getCountDays(Year) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Get the count of days in a given month of a given year of the Gregorian calendar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getYear()")
    call disp%show( getYear() )
    call disp%show("getMonth()")
    call disp%show( getMonth() )
    call disp%show("isLeapYear(getYear())")
    call disp%show( isLeapYear(getYear()) )
    call disp%show("getCountDays(getYear(), getMonth())")
    call disp%show( getCountDays(getYear(), getMonth()) )
    call disp%skip()

    call disp%skip()
    call disp%show("isLeapYear(2000_IK)")
    call disp%show( isLeapYear(2000_IK) )
    call disp%show("getCountDays(2000_IK, 2_IK)")
    call disp%show( getCountDays(2000_IK, 2_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("isLeapYear(2024_IK)")
    call disp%show( isLeapYear(2024_IK) )
    call disp%show("getCountDays(2024_IK, 2_IK)")
    call disp%show( getCountDays(2024_IK, 2_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("Year = [1_IK, 2_IK, 3_IK, 4_IK]")
                    Year = [1_IK, 2_IK, 3_IK, 4_IK]
    call disp%show("Year")
    call disp%show( Year )
    call disp%show("isLeapYear(Year)")
    call disp%show( isLeapYear(Year) )
    call disp%show("getCountDays(Year, 2_IK)")
    call disp%show( getCountDays(Year, 2_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("Year = -[1_IK, 2_IK, 3_IK, 4_IK]")
                    Year = -[1_IK, 2_IK, 3_IK, 4_IK]
    call disp%show("Year")
    call disp%show( Year )
    call disp%show("isLeapYear(Year)")
    call disp%show( isLeapYear(Year) )
    call disp%show("getCountDays(Year, 2_IK)")
    call disp%show( getCountDays(Year, 2_IK) )
    call disp%skip()

end program example