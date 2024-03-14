program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: ls, isFailedMakeDir

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("ls(SK_'.')")
    call disp%show( ls(SK_'.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("ls(SK_'..')")
    call disp%show( ls(SK_'..') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("ls(SK_'./*.F90')")
    call disp%show( ls(SK_'./*.F90') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("ls(SK_'./main.F90')")
    call disp%show( ls(SK_'./main.F90') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("ls(SK_'./*.nonexistent')")
    call disp%show( ls(SK_'./*.nonexistent') , deliml = SK_"""" )
    call disp%skip()

end program example