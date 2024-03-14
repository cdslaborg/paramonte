program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysShell, only: isFailedGetDirTemp

    implicit none

    logical(LK) :: failed
    character(255, SK) :: errmsg = SK_""
    character(:, SK), allocatable :: dirTemp

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("if (isFailedGetDirTemp(dirTemp)) then; call disp%show( 'Failed to fetch system shell temp. dir.' , deliml = SK_'''' ); else; call disp%show( dirTemp , deliml = SK_'''' ); end if")
                    if (isFailedGetDirTemp(dirTemp)) then; call disp%show( 'Failed to fetch system shell temp. dir.' , deliml = SK_'''' ); else; call disp%show( dirTemp , deliml = SK_'''' ); end if
    call disp%skip()

    call disp%skip()
    call disp%show("if (isFailedGetDirTemp(dirTemp, errmsg)) then; call disp%show( trim(errmsg) , deliml = SK_'''' ); else; call disp%show( dirTemp , deliml = SK_'''' ); end if")
                    if (isFailedGetDirTemp(dirTemp, errmsg)) then; call disp%show( trim(errmsg) , deliml = SK_'''' ); else; call disp%show( dirTemp , deliml = SK_'''' ); end if
    call disp%skip()

end program example