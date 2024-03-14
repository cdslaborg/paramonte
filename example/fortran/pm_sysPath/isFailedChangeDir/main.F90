program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: isFailedChangeDir, getDirCurrent, getDirSep

    implicit none

    character(:, SK), allocatable   :: cpath, npath
    logical(LK)                     :: failed

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("cpath = getDirCurrent() ! get the current working directory.")
                    cpath = getDirCurrent()
    call disp%show("cpath")
    call disp%show( cpath , deliml = SK_"""" )
    call disp%skip()
    call disp%show("npath = cpath//getDirSep()//SK_'..' ! new working directory.")
                    npath = cpath//getDirSep()//SK_'..'
    call disp%show("npath")
    call disp%show( npath , deliml = SK_"""" )
    call disp%skip()
    call disp%show("isFailedChangeDir(npath)")
    call disp%show( isFailedChangeDir(npath) )
    call disp%skip()
    call disp%show("getDirCurrent()")
    call disp%show( getDirCurrent() , deliml = SK_"""" )
    call disp%skip()

end program example