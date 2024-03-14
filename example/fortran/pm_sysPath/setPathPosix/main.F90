program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_sysPath, only: setPathPosix

    implicit none

    character(:), allocatable :: path

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("path = ''")
                    path = ''
    call disp%show("call setPathPosix(path)")
                    call setPathPosix(path)
    call disp%show("path")
    call disp%show( path , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'foo\bar'")
                    path = 'foo\bar'
    call disp%show("call setPathPosix(path)")
                    call setPathPosix(path)
    call disp%show("path")
    call disp%show( path , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'foo\/bar'")
                    path = 'foo\/bar'
    call disp%show("call setPathPosix(path)")
                    call setPathPosix(path)
    call disp%show("path")
    call disp%show( path , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'foo\/bar/'")
                    path = 'foo\/bar/'
    call disp%show("call setPathPosix(path)")
                    call setPathPosix(path)
    call disp%show("path")
    call disp%show( path , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'foo\/bar/\'")
                    path = 'foo\/bar/\'
    call disp%show("call setPathPosix(path)")
                    call setPathPosix(path)
    call disp%show("path")
    call disp%show( path , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'foo\/bar$/\'")
                    path = 'foo\/bar$/\'
    call disp%show("call setPathPosix(path)")
                    call setPathPosix(path)
    call disp%show("path")
    call disp%show( path , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = 'd:\foo\/bar$\/\'")
                    path = 'd:\foo\/bar$\/\'
    call disp%show("call setPathPosix(path, ignore = '\/')")
                    call setPathPosix(path, ignore = '\/')
    call disp%show("path")
    call disp%show( path , deliml = """" )
    call disp%skip()

end program example