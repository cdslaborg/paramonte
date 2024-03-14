program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_sysShell, only: isShellWindows
    use pm_sysShell, only: isShellPosix
    use pm_sysPath, only: getPathAbs

    implicit none

    character(:, SK), allocatable   :: path
    character(2047, SK)             :: errmsg
    logical(LK)                     :: failed

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isShellPosix()")
    call disp%show( isShellPosix() )
    call disp%show("isShellWindows()")
    call disp%show( isShellWindows() )
    call disp%skip()

    call disp%skip()
    call disp%show("path = ''")
                    path = ''
    call disp%show("path = getPathAbs(path)")
                    path = getPathAbs(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'foo\bar'")
                    path = 'foo\bar'
    call disp%show("path = getPathAbs(path)")
                    path = getPathAbs(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = '\foo\bar'")
                    path = '\foo\bar'
    call disp%show("path = getPathAbs(path)")
                    path = getPathAbs(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = '/foo\bar'")
                    path = '/foo\bar'
    call disp%show("path = getPathAbs(path)")
                    path = getPathAbs(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = './foo\/bar'")
                    path = './foo\/bar'
    call disp%show("path = getPathAbs(path)")
                    path = getPathAbs(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = '.\foo\/bar/'")
                    path = '.\foo\/bar/'
    call disp%show("path = getPathAbs(path)")
                    path = getPathAbs(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = '../foo\/bar/..'")
                    path = '../foo\/bar/..'
    call disp%show("path = getPathAbs(path)")
                    path = getPathAbs(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'foo\/bar/\'")
                    path = 'foo\/bar/\'
    call disp%show("path = getPathAbs(path)")
                    path = getPathAbs(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'd:\foo\/bar$\/\'")
                    path = 'd:\foo\/bar$\/\'
    call disp%show("path = getPathAbs(path)")
                    path = getPathAbs(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = '\\wsl$/Ubuntu-20.04' ! UNC path (recognized only by Windows OS).")
                    path = '\\wsl$/Ubuntu-20.04' ! UNC path (recognized only by Windows OS).
    call disp%show("path = getPathAbs(path)")
                    path = getPathAbs(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = '//wsl$/Ubuntu-20.04' ! UNC path (recognized by both Windows and POSIX OS).")
                    path = '//wsl$/Ubuntu-20.04' ! UNC path (recognized by both Windows and POSIX OS).
    call disp%show("path = getPathAbs(path)")
                    path = getPathAbs(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = getPathAbs(SK_"", failed, errmsg) ! handle potential runtime failures gracefully.")
                    path = getPathAbs(SK_"", failed, errmsg)
    call disp%show("failed")
    call disp%show( failed )
    call disp%skip()

end program example