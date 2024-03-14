program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysShell, only: isShellWindows

    implicit none

    character(255, SK)  :: errmsg = SK_""
    logical(LK)         :: failed

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isShellWindows()")
    call disp%show( isShellWindows() )
    call disp%skip()

    call disp%skip()
    call disp%show("isShellWindows(failed)")
    call disp%show( isShellWindows(failed) )
    call disp%show("failed")
    call disp%show( failed )
    call disp%skip()

    call disp%skip()
    call disp%show("isShellWindows(failed, errmsg)")
    call disp%show( isShellWindows(failed, errmsg) )
    call disp%show("failed")
    call disp%show( failed )
    call disp%show("trim(errmsg)")
    call disp%show( trim(errmsg) , deliml = SK_"""" )
    call disp%skip()

end program example