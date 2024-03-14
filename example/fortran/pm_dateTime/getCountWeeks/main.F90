program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: getCountWeeks
    use pm_arrayRange, only: getRange
    use pm_dateTime, only: getYear

    implicit none

    integer(IK), allocatable :: Year(:)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Get the count of weeks (with majority of their days) in the given year of the Gregorian calendar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getYear()")
    call disp%show( getYear() )
    call disp%show("getCountWeeks(getYear())")
    call disp%show( getCountWeeks(getYear()) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountWeeks(getRange(2000_IK, 2020_IK))")
    call disp%show( getCountWeeks(getRange(2000_IK, 2020_IK)) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Get the count of weeks (with majority of their days) in the given month of the year of the Gregorian calendar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getCountWeeks(getRange(2000_IK, 2020_IK), 1_IK)")
    call disp%show( getCountWeeks(getRange(2000_IK, 2020_IK), 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountWeeks(getRange(2000_IK, 2020_IK), 2_IK)")
    call disp%show( getCountWeeks(getRange(2000_IK, 2020_IK), 2_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountWeeks(getRange(2000_IK, 2020_IK), 3_IK)")
    call disp%show( getCountWeeks(getRange(2000_IK, 2020_IK), 3_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountWeeks(getRange(2000_IK, 2020_IK), 4_IK)")
    call disp%show( getCountWeeks(getRange(2000_IK, 2020_IK), 4_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountWeeks(getRange(2000_IK, 2020_IK), 5_IK)")
    call disp%show( getCountWeeks(getRange(2000_IK, 2020_IK), 5_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountWeeks(getRange(2000_IK, 2020_IK), 6_IK)")
    call disp%show( getCountWeeks(getRange(2000_IK, 2020_IK), 6_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountWeeks(getRange(2000_IK, 2020_IK), 7_IK)")
    call disp%show( getCountWeeks(getRange(2000_IK, 2020_IK), 7_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountWeeks(getRange(2000_IK, 2020_IK), 8_IK)")
    call disp%show( getCountWeeks(getRange(2000_IK, 2020_IK), 8_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountWeeks(getRange(2000_IK, 2020_IK), 8_IK)")
    call disp%show( getCountWeeks(getRange(2000_IK, 2020_IK), 8_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountWeeks(getRange(2000_IK, 2020_IK), 9_IK)")
    call disp%show( getCountWeeks(getRange(2000_IK, 2020_IK), 9_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountWeeks(getRange(2000_IK, 2020_IK), 10_IK)")
    call disp%show( getCountWeeks(getRange(2000_IK, 2020_IK), 10_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountWeeks(getRange(2000_IK, 2020_IK), 11_IK)")
    call disp%show( getCountWeeks(getRange(2000_IK, 2020_IK), 11_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountWeeks(getRange(2000_IK, 2020_IK), 12_IK)")
    call disp%show( getCountWeeks(getRange(2000_IK, 2020_IK), 12_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountWeeks(2022_IK, getRange(1_IK, 12_IK))")
    call disp%show( getCountWeeks(2022_IK, getRange(1_IK, 12_IK)) )
    call disp%skip()

end program example