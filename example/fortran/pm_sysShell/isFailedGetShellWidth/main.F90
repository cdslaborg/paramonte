program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysShell, only: isFailedGetShellWidth

    implicit none

    integer(IK) :: width
    logical(LK) :: failed

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show('if (.not. isFailedGetShellWidth(width)) call disp%show(width)')
                    if (.not. isFailedGetShellWidth(width)) call disp%show(width)
    call disp%skip()

end program example