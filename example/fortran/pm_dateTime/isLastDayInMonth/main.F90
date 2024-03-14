program example

    use pm_kind, only: SK, IK
    use pm_val2str, only: getStr
    use pm_io, only: display_type
    use pm_dateTime, only: isLastDayInMonth

    implicit none

    integer(IK) :: Values(8)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Determine if a given Gregorian calendar date is the last of the month.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("isLastDayInMonth(1582_IK, 10_IK, 15_IK)")
    call disp%show( isLastDayInMonth(1582_IK, 10_IK, 15_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("isLastDayInMonth(1_IK, 1_IK, 31_IK)")
    call disp%show( isLastDayInMonth(1_IK, 1_IK, 31_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("isLastDayInMonth(1999_IK, 12_IK, 31_IK)")
    call disp%show( isLastDayInMonth(1999_IK, 12_IK, 31_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("call date_and_time(values = Values(1:8))")
                    call date_and_time(values = Values(1:8))
    call disp%show("SK_'The statement `Today is the last day of month` is '//getStr(isLastDayInMonth(Values(1:3)))//SK_'.'")
    call disp%show( SK_'The statement `Today is the last day of month` is '//getStr(isLastDayInMonth(Values(1:3)))//SK_'.' , deliml = SK_"""" )
    call disp%skip()

end program example