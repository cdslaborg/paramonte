program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isCharUpper

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isCharUpper(SK_'A')")
    call disp%show( isCharUpper(SK_'A') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharUpper(SK_'Z')")
    call disp%show( isCharUpper(SK_'Z') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharUpper([character(1,SK) :: 'A', 'Z', 'L'])")
    call disp%show( isCharUpper([character(1,SK) :: 'A', 'Z', 'L']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharUpper([character(1,SK) :: 'a', 'Q', ' '])")
    call disp%show( isCharUpper([character(1,SK) :: 'a', 'Q', ' ']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharUpper(SK_' ')")
    call disp%show( isCharUpper(SK_' ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharUpper(SK_'1')")
    call disp%show( isCharUpper(SK_'1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharUpper(SK_'a')")
    call disp%show( isCharUpper(SK_'a') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharUpper(SK_'?')")
    call disp%show( isCharUpper(SK_'?') )
    call disp%skip()

end program example