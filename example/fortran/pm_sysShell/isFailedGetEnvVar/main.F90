program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysShell, only: isFailedGetEnvVar

    implicit none

    logical(LK) :: failed
    character(255, SK) :: errmsg = SK_""
    character(:, SK), allocatable :: value

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("if (isFailedGetEnvVar(SK_'OS', value)) then; call disp%show( 'Failed to fetch env. variable.' , deliml = SK_'''' ); else; call disp%show( value , deliml = SK_'''' ); end if")
                    if (isFailedGetEnvVar(SK_'OS', value)) then; call disp%show( 'Failed to fetch env. variable.' , deliml = SK_'''' ); else; call disp%show( value , deliml = SK_'''' ); end if
    call disp%skip()

    call disp%skip()
    call disp%show("if (isFailedGetEnvVar(SK_'PWD', value)) then; call disp%show( 'Failed to fetch env. variable.' , deliml = SK_'''' ); else; call disp%show( value , deliml = SK_'''' ); end if")
                    if (isFailedGetEnvVar(SK_'PWD', value)) then; call disp%show( 'Failed to fetch env. variable.' , deliml = SK_'''' ); else; call disp%show( value , deliml = SK_'''' ); end if
    call disp%skip()

    call disp%skip()
    call disp%show("if (isFailedGetEnvVar(SK_'HOME', value)) then; call disp%show( 'Failed to fetch env. variable.' , deliml = SK_'''' ); else; call disp%show( value , deliml = SK_'''' ); end if")
                    if (isFailedGetEnvVar(SK_'HOME', value)) then; call disp%show( 'Failed to fetch env. variable.' , deliml = SK_'''' ); else; call disp%show( value , deliml = SK_'''' ); end if
    call disp%skip()

    call disp%skip()
    call disp%show("if (isFailedGetEnvVar(SK_'SHELL', value)) then; call disp%show( 'Failed to fetch env. variable.' , deliml = SK_'''' ); else; call disp%show( value , deliml = SK_'''' ); end if")
                    if (isFailedGetEnvVar(SK_'SHELL', value)) then; call disp%show( 'Failed to fetch env. variable.' , deliml = SK_'''' ); else; call disp%show( value , deliml = SK_'''' ); end if
    call disp%skip()

    call disp%skip()
    call disp%show("if (isFailedGetEnvVar(SK_'', value, errmsg)) then; call disp%show( trim(errmsg) , deliml = SK_'''' ); else; call disp%show( value , deliml = SK_'''' ); end if")
                    if (isFailedGetEnvVar(SK_'', value, errmsg)) then; call disp%show( trim(errmsg) , deliml = SK_'''' ); else; call disp%show( value , deliml = SK_'''' ); end if
    call disp%skip()

end program example