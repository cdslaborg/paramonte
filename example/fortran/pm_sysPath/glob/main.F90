program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: glob, isFailedMakeDir

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("glob(SK_'.')")
    call disp%show( glob(SK_'.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("glob(SK_'../get*')")
    call disp%show( glob(SK_'../get*') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("glob(SK_'../get*/*')")
    call disp%show( glob(SK_'../get*/*') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("glob(SK_'./*.F90')")
    call disp%show( glob(SK_'./*.F90') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("glob(SK_'./*.nonexistent')")
    call disp%show( glob(SK_'./*.nonexistent') , deliml = SK_"""" )
    call disp%skip()

end program example