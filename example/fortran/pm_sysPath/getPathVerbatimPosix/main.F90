program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_sysPath, only: getPathVerbatimPosix

    implicit none

    character(:), allocatable :: path

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("path = ''")
                    path = ''
    call disp%show("path = getPathVerbatimPosix(path)")
                    path = getPathVerbatimPosix(path)
    call disp%show("path")
    call disp%show( path , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = './""paramonte""/library'")
                    path = "./""paramonte""/library"
    call disp%show("path = getPathVerbatimPosix(path)")
                    path = getPathVerbatimPosix(path)
    call disp%show("path")
    call disp%show( path , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show('path = "./''paramonte''/library"')
                    path = './''paramonte''/library'
    call disp%show("path = getPathVerbatimPosix(path)")
                    path = getPathVerbatimPosix(path)
    call disp%show("path")
    call disp%show( path , deliml = """" )
    call disp%skip()

end program example