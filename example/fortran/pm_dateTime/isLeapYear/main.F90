program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: isLeapYear, getYear

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getYear()")
    call disp%show( getYear() )
    call disp%show("isLeapYear(getYear())")
    call disp%show( isLeapYear(getYear()) )
    call disp%skip()

    call disp%skip()
    call disp%show("isLeapYear(2000_IK)")
    call disp%show( isLeapYear(2000_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("isLeapYear(2024_IK)")
    call disp%show( isLeapYear(2024_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("isLeapYear([1_IK, 2_IK, 3_IK, 4_IK])")
    call disp%show( isLeapYear([1_IK, 2_IK, 3_IK, 4_IK]) )
    call disp%skip()

    call disp%skip()
    call disp%show("isLeapYear(-[1_IK, 2_IK, 3_IK, 4_IK])")
    call disp%show( isLeapYear(-[1_IK, 2_IK, 3_IK, 4_IK]) )
    call disp%skip()

end program example