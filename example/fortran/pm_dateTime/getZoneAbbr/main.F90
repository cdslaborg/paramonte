program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: getZoneAbbr
    use pm_val2str, only: getStr

    implicit none

    integer(IK) :: zone
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getZoneAbbr()")
    call disp%show( getZoneAbbr() )
    call disp%skip()

    call disp%skip()
    do zone = -12 * 60 - 40, 14 * 60 + 40, 30
        call disp%skip()
        call disp%show("zone")
        call disp%show( zone )
        call disp%show("getZoneAbbr(zone)")
        call disp%show( getZoneAbbr(zone) , deliml = SK_"""" )
        call disp%skip()
    end do

end program example