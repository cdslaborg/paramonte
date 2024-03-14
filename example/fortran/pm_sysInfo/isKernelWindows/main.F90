program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysInfo, only: isKernelWindows

    implicit none

    logical(LK) :: failed
    type(display_type) :: disp
    character(255, SK) :: errmsg = SK_""

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isKernelWindows()")
    call disp%show( isKernelWindows() )
    call disp%skip()

    call disp%skip()
    call disp%show("isKernelWindows(failed)")
    call disp%show( isKernelWindows(failed) )
    call disp%show("failed")
    call disp%show( failed )
    call disp%skip()

    call disp%skip()
    call disp%show("isKernelWindows(failed, errmsg)")
    call disp%show( isKernelWindows(failed, errmsg) )
    call disp%show("failed")
    call disp%show( failed )
    call disp%show("trim(errmsg)")
    call disp%show( trim(errmsg) , deliml = SK_"""" )
    call disp%skip()

end program example