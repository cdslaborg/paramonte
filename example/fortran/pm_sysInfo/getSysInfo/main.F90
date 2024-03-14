program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysInfo, only: getSysInfo

    implicit none

    logical(LK) :: failed
    character(255, SK)  :: errmsg
    character(:, SK), allocatable :: sysinfo

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("sysinfo = getSysInfo()")
                    sysinfo = getSysInfo()
    call disp%show( sysinfo, deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("sysinfo = getSysInfo(failed)")
                    sysinfo = getSysInfo(failed)
    call disp%show( sysinfo, deliml = SK_"""" )
    call disp%show("failed ! Check if any error has occurred.")
    call disp%show( failed )
    if (failed) then
        call disp%show("SK_'error occurred.'")
        call disp%show( SK_'error occurred.' , deliml = SK_"""" )
    end if
    call disp%skip()

    call disp%skip()
    call disp%show("sysinfo = getSysInfo(failed, errmsg)")
                    sysinfo = getSysInfo(failed, errmsg)
    call disp%show( sysinfo, deliml = SK_"""" )
    call disp%show("failed ! Check if any error has occurred.")
    call disp%show( failed )
    if (failed) then
        call disp%show("SK_'error occurred: '//trim(errmsg)")
        call disp%show( SK_'error occurred: '//trim(errmsg) , deliml = SK_"""" )
    end if
    call disp%skip()

end program example