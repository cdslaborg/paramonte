program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: getZone
    use pm_val2str, only: getStr

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("SK_'The local time difference with respect to UTC is '//getStr(getZone())//SK_' minutes.'")
    call disp%show( SK_'The local time difference with respect to UTC is '//getStr(getZone())//SK_' minutes.' )
    call disp%skip()

end program example