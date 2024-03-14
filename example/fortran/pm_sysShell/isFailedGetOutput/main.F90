program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysShell, only: isFailedGetOutput

    implicit none

    logical(LK) :: failed
    character(255, SK) :: errmsg = SK_""
    character(:, SK), allocatable :: output

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("failed = isFailedGetOutput(SK_'ls', output) ! List all files in the current directory (works on POSIX-compliant systems).")
                    failed = isFailedGetOutput(SK_'ls', output)
    call disp%show("failed ! Check if any error has occurred.")
    call disp%show( failed )
    call disp%show("if (.not. failed) call disp%show( output , deliml = SK_'''' )")
                    if (.not. failed) call disp%show( output , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("failed = isFailedGetOutput(SK_'dir', output) ! List all files in the current directory (works on Windows-compliant systems).")
                    failed = isFailedGetOutput(SK_'dir', output)
    call disp%show("failed ! Check if any error has occurred.")
    call disp%show( failed )
    call disp%show("if (.not. failed) call disp%show( output , deliml = SK_'''' )")
                    if (.not. failed) call disp%show( output , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("failed = isFailedGetOutput(SK_'if shell syntax error happens, output remains unallocated.', output, errmsg).")
                    failed = isFailedGetOutput(SK_'if shell syntax error happens, output remains unallocated.', output, errmsg)
    call disp%show("failed ! Check if any error has occurred.")
    call disp%show( failed )
    call disp%show("if (failed) then; call disp%show( trim(errmsg) , deliml = SK_'''' ); else; call disp%show( output , deliml = SK_'''' ); end if")
                    if (failed) then; call disp%show( trim(errmsg) , deliml = SK_'''' ); else; call disp%show( output , deliml = SK_'''' ); end if
    call disp%skip()

    call disp%skip()
    call disp%show("failed = isFailedGetOutput(SK_'ls', output)")
                    failed = isFailedGetOutput(SK_'ls', output)
    call disp%show("failed ! Check if any error has occurred.")
    call disp%show( failed )
    call disp%show("if (.not. failed) call disp%show( output , deliml = SK_'''' )")
                    if (.not. failed) call disp%show( output , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("failed = isFailedGetOutput(SK_'ls -lhtr', output, errmsg) ! Gracefully handle runtime errors and capture the error message.")
                    failed = isFailedGetOutput(SK_'ls -lhtr', output, errmsg)
    call disp%show("failed ! Check if any error has occurred.")
    call disp%show( failed )
    call disp%show("if (failed) then; call disp%show( trim(errmsg) , deliml = SK_'''' ); else; call disp%show( output , deliml = SK_'''' ); end if")
                    if (failed) then; call disp%show( trim(errmsg) , deliml = SK_'''' ); else; call disp%show( output , deliml = SK_'''' ); end if
    call disp%skip()

end program example