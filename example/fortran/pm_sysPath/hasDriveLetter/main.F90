program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: hasDriveLetter

    implicit none

    character(1023, SK) :: errmsg = SK_""

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("hasDriveLetter('')")
    call disp%show( hasDriveLetter('') )
    call disp%skip()

    call disp%skip()
    call disp%show("hasDriveLetter('main.F90')")
    call disp%show( hasDriveLetter('main.F90') )
    call disp%skip()

    call disp%skip()
    call disp%show("hasDriveLetter('C:main.F90')")
    call disp%show( hasDriveLetter('C:main.F90') )
    call disp%skip()

    call disp%skip()
    call disp%show("hasDriveLetter(['/main.F90', '\main.F90'])")
    call disp%show( hasDriveLetter(['/main.F90', '\main.F90']) )
    call disp%skip()

    call disp%skip()
    call disp%show("hasDriveLetter(['C:/main.F90', 'C:\main.F90'])")
    call disp%show( hasDriveLetter(['C:/main.F90', 'C:\main.F90']) )
    call disp%skip()

    call disp%skip()
    call disp%show("hasDriveLetter(['../', '..\', '\..', '/..'])")
    call disp%show( hasDriveLetter(['../', '..\', '\..', '/..']) )
    call disp%skip()

    call disp%skip()
    call disp%show("hasDriveLetter(['C:\', '..\', '\..', '/..'])")
    call disp%show( hasDriveLetter(['C:\', '..\', '\..', '/..']) )
    call disp%skip()

    call disp%skip()
    call disp%show("hasDriveLetter(['.', '\', '/'])")
    call disp%show( hasDriveLetter(['.', '\', '/']) )
    call disp%skip()

end program example