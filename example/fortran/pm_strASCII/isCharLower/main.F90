program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isCharLower

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isCharLower(SK_'a')")
    call disp%show( isCharLower(SK_'a') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharLower(SK_'z')")
    call disp%show( isCharLower(SK_'z') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharLower([character(1,SK) :: 'a', 'z', 'l'])")
    call disp%show( isCharLower([character(1,SK) :: 'a', 'z', 'l']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharLower([character(1,SK) :: 'A', 'q', ' '])")
    call disp%show( isCharLower([character(1,SK) :: 'A', 'q', ' ']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharLower(SK_' ')")
    call disp%show( isCharLower(SK_' ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharLower(SK_'1')")
    call disp%show( isCharLower(SK_'1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharLower(SK_'A')")
    call disp%show( isCharLower(SK_'A') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharLower(SK_'?')")
    call disp%show( isCharLower(SK_'?') )
    call disp%skip()

end program example