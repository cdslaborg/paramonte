program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isStrAlphaNumAll

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isStrAlphaNumAll('1')")
    call disp%show( isStrAlphaNumAll('1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAll('0123')")
    call disp%show( isStrAlphaNumAll('0123') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAll('a01C3F')")
    call disp%show( isStrAlphaNumAll('a01C3F') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAll('ParaMonte')")
    call disp%show( isStrAlphaNumAll('ParaMonte') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAll([character(4) :: '1a2b3c', 'paramonte', '0123'])")
    call disp%show( isStrAlphaNumAll([character(4) :: '1a2b3c', 'paramonte', '0123']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAll([character(4) :: ' ', '1aBc2', '#@!?'])")
    call disp%show( isStrAlphaNumAll([character(4) :: ' ', '1aBc2', '#@!?']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAll('-1')")
    call disp%show( isStrAlphaNumAll('-1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAll('+1')")
    call disp%show( isStrAlphaNumAll('+1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAll('1.')")
    call disp%show( isStrAlphaNumAll('1.') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAll('')")
    call disp%show( isStrAlphaNumAll('') )
    call disp%skip()

end program example