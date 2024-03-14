program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: isFile

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isFile(SK_'')")
    call disp%show( isFile(SK_'') )
    call disp%skip()

    call disp%skip()
    call disp%show("isFile(SK_'main.F90')")
    call disp%show( isFile(SK_'main.F90') )
    call disp%skip()

    call disp%skip()
    call disp%show("isFile(SK_'main.out.F90')")
    call disp%show( isFile(SK_'main.out.F90') )
    call disp%skip()

    call disp%skip()
    call disp%show("isFile(SK_'paramonte')")
    call disp%show( isFile(SK_'paramonte') )
    call disp%skip()

    call disp%skip()
    call disp%show("isFile(SK_'.')")
    call disp%show( isFile(SK_'.') )
    call disp%skip()

    call disp%skip()
    call disp%show("isFile(SK_'./')")
    call disp%show( isFile(SK_'./') )
    call disp%skip()

    call disp%skip()
    call disp%show("isFile(SK_'..')")
    call disp%show( isFile(SK_'..') )
    call disp%skip()

end program example