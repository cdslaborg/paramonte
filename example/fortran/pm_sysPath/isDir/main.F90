program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: isDir, isFailedMakeDir

    implicit none

    character(1023, SK) :: errmsg = SK_""

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isDir(SK_'')")
    call disp%show( isDir(SK_'') )
    call disp%skip()

    call disp%skip()
    call disp%show("isDir(SK_'main.F90')")
    call disp%show( isDir(SK_'main.F90') )
    call disp%skip()

    call disp%skip()
    call disp%show("isDir(SK_'./paramonte')")
    call disp%show( isDir(SK_'./paramonte') )
    call disp%show("isFailedMakeDir(SK_'./paramonte')")
    call disp%show( isFailedMakeDir(SK_'./paramonte') )
    call disp%show("isDir(SK_'./paramonte/')")
    call disp%show( isDir(SK_'./paramonte/') )
    call disp%show("isDir(SK_'./paramonte')")
    call disp%show( isDir(SK_'./paramonte') )
    call disp%skip()

    call disp%skip()
    call disp%show("isDir(SK_'.')")
    call disp%show( isDir(SK_'.') )
    call disp%skip()

    call disp%skip()
    call disp%show("isDir(SK_'./')")
    call disp%show( isDir(SK_'./') )
    call disp%skip()

    call disp%skip()
    call disp%show("isDir(SK_'..')")
    call disp%show( isDir(SK_'..') )
    call disp%skip()

end program example