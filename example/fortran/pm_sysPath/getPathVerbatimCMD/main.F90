program example

    use pm_kind, only: SK, LK
    use pm_io, only: display_type
    use pm_sysPath, only: getPathVerbatimCMD

    implicit none

    character(:), allocatable :: path

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("path = ''")
                    path = ''
    call disp%show("path = getPathVerbatimCMD(path)")
                    path = getPathVerbatimCMD(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = './paramonte'")
                    path = './paramonte'
    call disp%show("path = getPathVerbatimCMD(path)")
                    path = getPathVerbatimCMD(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("path = './paramonte""'")
                    path = './paramonte"'
    call disp%show("path = getPathVerbatimCMD(path)")
                    path = getPathVerbatimCMD(path)
    call disp%show("path")
    call disp%show( path , deliml = SK_'''' )
    call disp%skip()

end program example