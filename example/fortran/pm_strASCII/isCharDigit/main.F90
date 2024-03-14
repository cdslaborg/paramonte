program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isCharDigit

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isCharDigit(SK_'1')")
    call disp%show( isCharDigit(SK_'1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharDigit(SK_'-')")
    call disp%show( isCharDigit(SK_'-') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharDigit(SK_' ')")
    call disp%show( isCharDigit(SK_' ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharDigit([character(1,SK) :: '1', 'p', '0'])")
    call disp%show( isCharDigit([character(1,SK) :: '1', 'p', '0']) )
    call disp%skip()

end program example