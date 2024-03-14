program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysShell, only: shellis_type

    implicit none

    type(shellis_type)  :: shellis
    logical(LK)         :: failed
    character(255, SK)  :: errmsg

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("shellis = shellis_type()")
                    shellis = shellis_type()
    call dispShellIs()
    call disp%skip()

    call disp%skip()
    call disp%show("shellis = shellis_type(failed)")
                    shellis = shellis_type(failed)
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
    call disp%show("shellis = shellis_type(failed, errmsg)")
                    shellis = shellis_type(failed, errmsg)
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
        call disp%show("shellis%ash         ")
        call disp%show( shellis%ash          )
        call disp%show("shellis%bash        ")
        call disp%show( shellis%bash         )
        call disp%show("shellis%csh         ")
        call disp%show( shellis%csh          )
        call disp%show("shellis%cmd         ")
        call disp%show( shellis%cmd          )
        call disp%show("shellis%dash        ")
        call disp%show( shellis%dash         )
        call disp%show("shellis%fish        ")
        call disp%show( shellis%fish         )
        call disp%show("shellis%ksh         ")
        call disp%show( shellis%ksh          )
        call disp%show("shellis%posix       ")
        call disp%show( shellis%posix        )
        call disp%show("shellis%powershell  ")
        call disp%show( shellis%powershell   )
        call disp%show("shellis%sh          ")
        call disp%show( shellis%sh           )
        call disp%show("shellis%tcsh        ")
        call disp%show( shellis%tcsh         )
        call disp%show("shellis%windows     ")
        call disp%show( shellis%windows      )
        call disp%show("shellis%zsh         ")
        call disp%show( shellis%zsh          )
        call disp%show("shellis%yash        ")
        call disp%show( shellis%yash         )
    end subroutine

end program example