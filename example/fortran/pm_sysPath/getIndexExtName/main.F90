program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_sysPath, only: getIndexExtName

    implicit none

    character(:), allocatable :: path

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("path = SK_''")
                    path = SK_''
    call disp%show("getIndexExtName(path, SK_'/')")
    call disp%show( getIndexExtName(path, SK_'/') )
    call disp%show("path(getIndexExtName(path, SK_'/'):)")
    call disp%show( path(getIndexExtName(path, SK_'/'):) , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'.'")
                    path = SK_'.'
    call disp%show("getIndexExtName(path, SK_'/')")
    call disp%show( getIndexExtName(path, SK_'/') )
    call disp%show("path(getIndexExtName(path, SK_'/'):)")
    call disp%show( path(getIndexExtName(path, SK_'/'):) , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'..'")
                    path = SK_'..'
    call disp%show("getIndexExtName(path, SK_'/')")
    call disp%show( getIndexExtName(path, SK_'/') )
    call disp%show("path(getIndexExtName(path, SK_'/'):)")
    call disp%show( path(getIndexExtName(path, SK_'/'):) , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'/..'")
                    path = SK_'/..'
    call disp%show("getIndexExtName(path, SK_'/')")
    call disp%show( getIndexExtName(path, SK_'/') )
    call disp%show("path(getIndexExtName(path, SK_'/'):)")
    call disp%show( path(getIndexExtName(path, SK_'/'):) , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'../'")
                    path = SK_'../'
    call disp%show("getIndexExtName(path, SK_'/')")
    call disp%show( getIndexExtName(path, SK_'/') )
    call disp%show("path(getIndexExtName(path, SK_'/'):)")
    call disp%show( path(getIndexExtName(path, SK_'/'):) , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'/'")
                    path = SK_'/'
    call disp%show("getIndexExtName(path, SK_'/')")
    call disp%show( getIndexExtName(path, SK_'/') )
    call disp%show("path(getIndexExtName(path, SK_'/'):)")
    call disp%show( path(getIndexExtName(path, SK_'/'):) , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'/paramonte'")
                    path = SK_'/paramonte'
    call disp%show("getIndexExtName(path, SK_'/')")
    call disp%show( getIndexExtName(path, SK_'/') )
    call disp%show("path(getIndexExtName(path, SK_'/'):)")
    call disp%show( path(getIndexExtName(path, SK_'/'):) , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'./paramonte'")
                    path = SK_'./paramonte'
    call disp%show("getIndexExtName(path, SK_'/')")
    call disp%show( getIndexExtName(path, SK_'/') )
    call disp%show("path(getIndexExtName(path, SK_'/'):)")
    call disp%show( path(getIndexExtName(path, SK_'/'):) , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'./paramonte/parallel/library.txt'")
                    path = SK_'./paramonte/parallel/library.txt'
    call disp%show("getIndexExtName(path, SK_'/')")
    call disp%show( getIndexExtName(path, SK_'/') )
    call disp%show("path(getIndexExtName(path, SK_'/'):)")
    call disp%show( path(getIndexExtName(path, SK_'/'):) , deliml = """" )
    call disp%skip()

end program example