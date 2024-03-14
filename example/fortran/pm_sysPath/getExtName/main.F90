program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_sysPath, only: getExtName

    implicit none

    character(:), allocatable :: path

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("path = SK_''")
                    path = SK_''
    call disp%show("getExtName(path, SK_'/')")
    call disp%show( getExtName(path, SK_'/') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'.'")
                    path = SK_'.'
    call disp%show("getExtName(path, SK_'/')")
    call disp%show( getExtName(path, SK_'/') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'..'")
                    path = SK_'..'
    call disp%show("getExtName(path, SK_'/')")
    call disp%show( getExtName(path, SK_'/') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'/..'")
                    path = SK_'/..'
    call disp%show("getExtName(path, SK_'/')")
    call disp%show( getExtName(path, SK_'/') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'../'")
                    path = SK_'../'
    call disp%show("getExtName(path, SK_'/')")
    call disp%show( getExtName(path, SK_'/') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'/'")
                    path = SK_'/'
    call disp%show("getExtName(path, SK_'/')")
    call disp%show( getExtName(path, SK_'/') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'/paramonte'")
                    path = SK_'/paramonte'
    call disp%show("getExtName(path, SK_'/')")
    call disp%show( getExtName(path, SK_'/') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'./paramonte'")
                    path = SK_'./paramonte'
    call disp%show("getExtName(path, SK_'/')")
    call disp%show( getExtName(path, SK_'/') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'./paramonte/parallel/library.txt'")
                    path = SK_'./paramonte/parallel/library.txt'
    call disp%show("getExtName(path, SK_'/')")
    call disp%show( getExtName(path, SK_'/') , deliml = SK_"""" )
    call disp%skip()

end program example