program example

    use pm_kind, only: SK, IK
    use pm_val2str, only: getStr
    use pm_io, only: display_type
    use pm_dateTime, only: getDateBefore

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
    call disp%show("getDateBefore(1582_IK, 10_IK, 15_IK)")
    call disp%show( getDateBefore(1582_IK, 10_IK, 15_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateBefore(1_IK, 1_IK, 1_IK)")
    call disp%show( getDateBefore(1_IK, 1_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateBefore(-1_IK, 1_IK, 1_IK)")
    call disp%show( getDateBefore(-1_IK, 1_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateBefore(1999_IK, 12_IK, 1_IK)")
    call disp%show( getDateBefore(1999_IK, 12_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateBefore(2000_IK, 1_IK, 1_IK)")
    call disp%show( getDateBefore(2000_IK, 1_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateBefore(2000_IK, 3_IK, 1_IK) ! leap year")
    call disp%show( getDateBefore(2000_IK, 3_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateBefore(1999_IK, 3_IK, 1_IK) ! non-leap year")
    call disp%show( getDateBefore(1999_IK, 3_IK, 1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("call date_and_time(values = Values(1:8))")
                    call date_and_time(values = Values(1:8))
    call disp%show("SK_'The statement `The Gregorian date before today ('//getStr(Values(1:3))//SK_') was ('//getStr(getDateBefore(Values(1:3)))//SK_').'")
    call disp%show( SK_'The statement `The Gregorian date before today ('//getStr(Values(1:3))//SK_') was ('//getStr(getDateBefore(Values(1:3)))//SK_').' , deliml = SK_"""" )
    call disp%skip()

end program example