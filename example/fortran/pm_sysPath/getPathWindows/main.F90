program example

    use pm_kind, only: SK, LK
    use pm_io, only: display_type
    use pm_sysPath, only: getPathWindows

    implicit none

    character(:), allocatable :: path

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("path = ''")
                    path = ''
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = '\ '")
                    path = '\ '
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = '\\\'")
                    path = '\\\'
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = '\\\./\'")
                    path = '\\\./\'
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = '/ '")
                    path = '/ '
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = '/ /'")
                    path = '/ /'
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = './Para Monte/'")
                    path = './Para Monte/'
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'foo\$bar'")
                    path = 'foo\$bar'
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'foo\\$bar'")
                    path = 'foo\\$bar'
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'foo/\$bar'")
                    path = 'foo/\$bar'
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'foo/\$bar '")
                    path = 'foo/\$bar '
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'foo\{/bar/'")
                    path = 'foo\{/bar/'
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'foo\/bar$/\'")
                    path = 'foo\/bar$/\'
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'd:/\foo\/bar$/\'")
                    path = 'd:/\foo\/bar$/\'
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'd:\\foo\/bar$/\'")
                    path = 'd:\\foo\/bar$/\'
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = './paramonte library'")
                    path = './paramonte library'
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = '\\wsl$\Ubuntu-20.04\'")
                    path = '\\wsl$\Ubuntu-20.04\'
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = '//wsl$/Ubuntu-20.04//'")
                    path = '//wsl$/Ubuntu-20.04//'
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = '""./paramonte library""'")
                    path = '"./paramonte library"'
    call disp%show("path = getPathWindows(path)")
                    path = getPathWindows(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'd:/foo\/bar$\//'")
                    path = 'd:/foo\/bar$\//'
    call disp%show("path = getPathWindows(path, ignore = '\/')")
                    path = getPathWindows(path, ignore = '\/')
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

end program example