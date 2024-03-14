program example

    use pm_kind, only: SK, LK
    use pm_io, only: display_type
    use pm_sysPath, only: getPathVerbatim
    use pm_sysShell, only: shell_type

    implicit none

    type(shell_type) :: Shell
    character(:, SK), allocatable :: path

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show('shell = shell_type()')
                    shell = shell_type()
    call disp%show("shell%name")
    call disp%show( shell%name , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_''")
                    path = SK_''
    call disp%show("path = getPathVerbatim(path)")
                    path = getPathVerbatim(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = './""paramonte""/library'")
                    path = "./""paramonte""/library"
    call disp%show("path = getPathVerbatim(path)")
                    path = getPathVerbatim(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show('path = "./''paramonte''/library"')
                    path = './''paramonte''/library'
    call disp%show("path = getPathVerbatim(path)")
                    path = getPathVerbatim(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show('path = "./''paramonte''/library"')
                    path = './''paramonte''/library'
    call disp%show("path = getPathVerbatim(path)")
                    path = getPathVerbatim(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_"""" )
    call disp%skip()

end program example