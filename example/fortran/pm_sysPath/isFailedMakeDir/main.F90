program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: isFailedMakeDir, isDir, getPathNew, getDirSep

    implicit none

    character(:, SK), allocatable :: cpath

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("cpath = getPathNew() ! define a random unique path name.")
                    cpath = getPathNew()
    call disp%show("cpath")
    call disp%show( cpath , deliml = SK_"""" )
    call disp%show("isDir(cpath)")
    call disp%show( isDir(cpath) )
    call disp%show("isFailedMakeDir(cpath)")
    call disp%show( isFailedMakeDir(cpath) )
    call disp%show("isDir(cpath)")
    call disp%show( isDir(cpath) )
    call disp%skip()

    call disp%skip()
    call disp%show("cpath = getPathNew(dir = getPathNew())) ! define a random unique **nested** path name.")
                    cpath = getPathNew(dir = getPathNew())
    call disp%show("cpath")
    call disp%show( cpath , deliml = SK_"""" )
    call disp%show("isDir(cpath)")
    call disp%show( isDir(cpath) )
    call disp%show("isFailedMakeDir(cpath)")
    call disp%show( isFailedMakeDir(cpath) )
    call disp%show("isDir(cpath)")
    call disp%show( isDir(cpath) )
    call disp%skip()

end program example