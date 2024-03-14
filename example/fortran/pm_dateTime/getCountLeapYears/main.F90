program example

    use pm_kind, only: SK, IK
    use pm_val2str, only: getStr
    use pm_io, only: display_type
    use pm_dateTime, only: getCountLeapYears
    use pm_arrayRange, only: getRange
    use pm_dateTime, only: getYear

    implicit none

    integer(IK), allocatable :: Year(:)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Determine if a given Gregorian calendar date is the last of the month.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getCountLeapYears(until = 1_IK)")
    call disp%show( getCountLeapYears(until = 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountLeapYears(until = 2000_IK)")
    call disp%show( getCountLeapYears(until = 2000_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getYear()")
    call disp%show( getYear() )
    call disp%show("getCountLeapYears(until = getYear())")
    call disp%show( getCountLeapYears(until = getYear()) )
    call disp%skip()

    call disp%skip()
    call disp%show("getYear()")
    call disp%show( getYear() )
    call disp%show("getCountLeapYears(until = getYear(), since = 1_IK)")
    call disp%show( getCountLeapYears(until = getYear(), since = 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountLeapYears(until = [1_IK, 2_IK, 3_IK, 4_IK])")
    call disp%show( getCountLeapYears(until = [1_IK, 2_IK, 3_IK, 4_IK]) )
    call disp%skip()

    call disp%skip()
    call disp%show("Year = getRange(-5_IK, 5_IK, 1_IK)")
                    Year = getRange(-5_IK, 5_IK, 1_IK)
    call disp%show("Year")
    call disp%show( Year )
    call disp%show("getCountLeapYears(until = Year, since = -5_IK)")
    call disp%show( getCountLeapYears(until = Year, since = -5_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("Year = getRange(-8_IK, 8_IK, 1_IK)")
                    Year = getRange(-8_IK, 8_IK, 1_IK)
    call disp%show("Year")
    call disp%show( Year )
    call disp%show("getCountLeapYears(until = Year)")
    call disp%show( getCountLeapYears(until = Year) )
    call disp%skip()

    call disp%skip()
    call disp%show("getYear()")
    call disp%show( getYear() )
    call disp%show("getCountLeapYears(until = getYear(), since = 1_IK)")
    call disp%show( getCountLeapYears(until = getYear(), since = 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getYear()")
    call disp%show( getYear() )
    call disp%show("getCountLeapYears(until = getYear(), since = [1999_IK, 2000_IK, 2001_IK])")
    call disp%show( getCountLeapYears(until = getYear(), since = [1999_IK, 2000_IK, 2001_IK]) )
    call disp%skip()

end program example