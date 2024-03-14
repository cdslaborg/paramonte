program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysInfo, only: kernel_type

    implicit none

    logical(LK) :: failed
    character(255, SK) :: errmsg
    type(kernel_type) :: kernel

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("kernel = kernel_type()")
                    kernel = kernel_type()
    call dispkernel()
    call disp%skip()

    call disp%skip()
    call disp%show("kernel = kernel_type(failed)")
                    kernel = kernel_type(failed)
    call disp%show("failed ! Check if any error has occurred.")
    call disp%show( failed )
    if (failed) then
    call disp%show("SK_'error occurred.'")
    call disp%show( SK_'error occurred.' , deliml = SK_"""" )
    else
    call dispkernel()
    end if
    call disp%skip()

    call disp%skip()
    call disp%show("kernel = kernel_type(failed, errmsg)")
                    kernel = kernel_type(failed, errmsg)
    call disp%show("failed ! Check if any error has occurred.")
    call disp%show( failed )
    if (failed) then
    call disp%show("errmsg")
    call disp%show( errmsg , deliml = SK_"""" )
    else
    call dispkernel()
    end if
    call disp%skip()

contains

    subroutine dispkernel()
#if     __INTEL_COMPILER
#undef  linux
#endif
        call disp%show("kernel%name")
        call disp%show( kernel%name , deliml = SK_"""" )
        call disp%show("kernel%is%windows")
        call disp%show( kernel%is%windows )
        call disp%show("kernel%is%cygwin")
        call disp%show( kernel%is%cygwin )
        call disp%show("kernel%is%mingw")
        call disp%show( kernel%is%mingw )
        call disp%show("kernel%is%msys")
        call disp%show( kernel%is%msys )
        call disp%show("kernel%is%linux")
        call disp%show( kernel%is%linux )
        call disp%show("kernel%is%darwin")
        call disp%show( kernel%is%darwin )
        call disp%show("kernel%is%freebsd")
        call disp%show( kernel%is%freebsd )
    end subroutine

end program example