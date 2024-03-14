program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysShell, only: shell_type

    implicit none

    type(shell_type)    :: Shell
    logical(LK)         :: failed
    character(255, SK)  :: errmsg

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("shell = shell_type()")
                    shell = shell_type()
    call dispShell()
    call disp%skip()

    call disp%skip()
    call disp%show("shell = shell_type(failed)")
                    shell = shell_type(failed)
    call disp%show("failed ! Check if any error has occurred.")
    call disp%show( failed )
    if (failed) then
    call disp%show("SK_'error occurred.'")
    call disp%show( SK_'error occurred.' , deliml = SK_"""" )
    else
    call dispShell()
    end if
    call disp%skip()

    call disp%skip()
    call disp%show("shell = shell_type(failed, errmsg)")
                    shell = shell_type(failed, errmsg)
    call disp%show("failed ! Check if any error has occurred.")
    call disp%show( failed )
    if (failed) then
    call disp%show("errmsg")
    call disp%show( errmsg , deliml = SK_"""" )
    else
    call dispShell()
    end if
    call disp%skip()

contains

    subroutine dispShell()
        call disp%show("shell%is%ash        ")
        call disp%show( shell%is%ash         )
        call disp%show("shell%is%bash       ")
        call disp%show( shell%is%bash        )
        call disp%show("shell%is%csh        ")
        call disp%show( shell%is%csh         )
        call disp%show("shell%is%cmd        ")
        call disp%show( shell%is%cmd         )
        call disp%show("shell%is%dash       ")
        call disp%show( shell%is%dash        )
        call disp%show("shell%is%fish       ")
        call disp%show( shell%is%fish        )
        call disp%show("shell%is%ksh        ")
        call disp%show( shell%is%ksh         )
        call disp%show("shell%is%posix      ")
        call disp%show( shell%is%posix       )
        call disp%show("shell%is%powershell ")
        call disp%show( shell%is%powershell  )
        call disp%show("shell%is%sh         ")
        call disp%show( shell%is%sh          )
        call disp%show("shell%is%tcsh       ")
        call disp%show( shell%is%tcsh        )
        call disp%show("shell%is%windows    ")
        call disp%show( shell%is%windows     )
        call disp%show("shell%is%zsh        ")
        call disp%show( shell%is%zsh         )
        call disp%show("shell%is%yash       ")
        call disp%show( shell%is%yash        )
        call disp%show("shell%dirsep        ")
        call disp%show( shell%dirsep         , deliml = SK_"""" )
        call disp%show("shell%name          ")
        call disp%show( shell%name           , deliml = SK_"""" )
    end subroutine

end program example