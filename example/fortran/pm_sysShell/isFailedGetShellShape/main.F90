program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysShell, only: isFailedGetShellShape

    implicit none

    logical(LK) :: failed
    integer(IK) :: width, height

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show('if (.not. isFailedGetShellShape(width, height)) call disp%show([width, height])')
                    if (.not. isFailedGetShellShape(width, height)) call disp%show([width, height])
    call disp%skip()

end program example