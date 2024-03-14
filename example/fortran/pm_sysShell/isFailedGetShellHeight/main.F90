program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysShell, only: isFailedGetShellHeight

    implicit none

    logical(LK) :: failed
    integer(IK) :: height

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show('if (.not. isFailedGetShellHeight(height)) call disp%show(height)')
                    if (.not. isFailedGetShellHeight(height)) call disp%show(height)
    call disp%skip()

end program example