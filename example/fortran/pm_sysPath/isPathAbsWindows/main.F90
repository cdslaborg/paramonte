program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: isPathAbsWindows

    implicit none

    character(1023, SK) :: errmsg = SK_""

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isPathAbsWindows('')")
    call disp%show( isPathAbsWindows('') )
    call disp%skip()

    call disp%skip()
    call disp%show("isPathAbsWindows(['.', '\', '/'])")
    call disp%show( isPathAbsWindows(['.', '\', '/']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isPathAbsWindows(['./', '\.', '/.'])")
    call disp%show( isPathAbsWindows(['./', '\.', '/.']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isPathAbsWindows(['../', '..\', '\..', '/..'])")
    call disp%show( isPathAbsWindows(['../', '..\', '\..', '/..']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isPathAbsWindows('main.F90')")
    call disp%show( isPathAbsWindows('main.F90') )
    call disp%skip()

    call disp%skip()
    call disp%show("isPathAbsWindows(['/main.F90', '\main.F90'])")
    call disp%show( isPathAbsWindows(['/main.F90', '\main.F90']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isPathAbsWindows('C:main.F90')")
    call disp%show( isPathAbsWindows('C:main.F90') )
    call disp%skip()

    call disp%skip()
    call disp%show("isPathAbsWindows(['C:/main.F90', 'C:\main.F90'])")
    call disp%show( isPathAbsWindows(['C:/main.F90', 'C:\main.F90']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isPathAbsWindows(['./paramonte', '//paramonte', '\/paramonte', '\\paramonte'])")
    call disp%show( isPathAbsWindows(['./paramonte', '//paramonte', '\/paramonte', '\\paramonte']) )
    call disp%skip()

end program example