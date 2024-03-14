program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysInfo, only: kernelis_type

    implicit none

    logical(LK) :: failed
    character(255, SK) :: errmsg
    type(kernelis_type) :: kernelis

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("kernelis = kernelis_type()")
                    kernelis = kernelis_type()
    call dispShellIs()
    call disp%skip()

    call disp%skip()
    call disp%show("kernelis = kernelis_type(failed)")
                    kernelis = kernelis_type(failed)
    call disp%show("failed ! Check if any error has occurred.")
    call disp%show( failed )
    if (failed) then
    call disp%show("SK_'error occurred.'")
    call disp%show( SK_'error occurred.' , deliml = SK_"""" )
    else
    call dispShellIs()
    end if
    call disp%skip()

    call disp%skip()
    call disp%show("kernelis = kernelis_type(failed, errmsg)")
                    kernelis = kernelis_type(failed, errmsg)
    call disp%show("failed ! Check if any error has occurred.")
    call disp%show( failed )
    if (failed) then
    call disp%show("errmsg")
    call disp%show( errmsg , deliml = SK_"""" )
    else
    call dispShellIs()
    end if
    call disp%skip()

contains

    subroutine dispShellIs()
#if     __INTEL_COMPILER
#undef  linux
#endif
        call disp%show("kernelis%windows")
        call disp%show( kernelis%windows )
        call disp%show("kernelis%cygwin")
        call disp%show( kernelis%cygwin )
        call disp%show("kernelis%mingw")
        call disp%show( kernelis%mingw )
        call disp%show("kernelis%msys")
        call disp%show( kernelis%msys )
        call disp%show("kernelis%linux")
        call disp%show( kernelis%linux )
        call disp%show("kernelis%darwin")
        call disp%show( kernelis%darwin )
        call disp%show("kernelis%freebsd")
        call disp%show( kernelis%freebsd )
    end subroutine

end program example