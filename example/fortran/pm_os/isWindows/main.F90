program example

    use pm_kind, only: LK, IK, SK
    use pm_os, only: isDarwin, isLinux, isWindows
    use pm_io, only: display_type
    use pm_val2int, only: getInt

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("[isDarwin(), isLinux(), isWindows()]")
    call disp%show( [isDarwin(), isLinux(), isWindows()] )
    call disp%skip()

    call disp%skip()
    call disp%show("if (sum(getInt([isDarwin(), isLinux(), isWindows()])) /= 1) error stop 'Multiple Operating Systems cannot exist simultaneously at runtime!'")
                    if (sum(getInt([isDarwin(), isLinux(), isWindows()])) /= 1) error stop 'Multiple Operating Systems cannot exist simultaneously at runtime!'
    call disp%skip()

end program example