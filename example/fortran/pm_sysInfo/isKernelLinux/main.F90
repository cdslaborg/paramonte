program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysInfo, only: isKernelLinux

    implicit none

    logical(LK) :: failed
    type(display_type) :: disp
    character(255, SK) :: errmsg = SK_""

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isKernelLinux()")
    call disp%show( isKernelLinux() )
    call disp%skip()

    call disp%skip()
    call disp%show("isKernelLinux(failed)")
    call disp%show( isKernelLinux(failed) )
    call disp%show("failed")
    call disp%show( failed )
    call disp%skip()

    call disp%skip()
    call disp%show("isKernelLinux(failed, errmsg)")
    call disp%show( isKernelLinux(failed, errmsg) )
    call disp%show("failed")
    call disp%show( failed )
    call disp%show("trim(errmsg)")
    call disp%show( trim(errmsg) , deliml = SK_"""" )
    call disp%skip()

end program example