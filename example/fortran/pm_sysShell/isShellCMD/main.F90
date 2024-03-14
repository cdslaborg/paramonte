program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysShell, only: isShellCMD

    implicit none

    character(255, SK)  :: errmsg = SK_""
    logical(LK)         :: failed

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isShellCMD()")
    call disp%show( isShellCMD() )
    call disp%skip()

    call disp%skip()
    call disp%show("isShellCMD(failed)")
    call disp%show( isShellCMD(failed) )
    call disp%show("failed")
    call disp%show( failed )
    call disp%skip()

    call disp%skip()
    call disp%show("isShellCMD(failed, errmsg)")
    call disp%show( isShellCMD(failed, errmsg) )
    call disp%show("failed")
    call disp%show( failed )
    call disp%show("trim(errmsg)")
    call disp%show( trim(errmsg) , deliml = SK_"""" )
    call disp%skip()

end program example