program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isStrAlphaNumAny

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isStrAlphaNumAny(SK_'1')")
    call disp%show( isStrAlphaNumAny(SK_'1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAny(SK_'0123')")
    call disp%show( isStrAlphaNumAny(SK_'0123') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAny(SK_'a01C3F')")
    call disp%show( isStrAlphaNumAny(SK_'a01C3F') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAny(SK_'ParaMonte')")
    call disp%show( isStrAlphaNumAny(SK_'ParaMonte') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAny([character(4,SK) :: '1a2b3c', 'paramonte', '0123'])")
    call disp%show( isStrAlphaNumAny([character(4,SK) :: '1a2b3c', 'paramonte', '0123']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAny([character(4,SK) :: ' ', '1aBc2', '#@!?'])")
    call disp%show( isStrAlphaNumAny([character(4,SK) :: ' ', '1aBc2', '#@!?']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAny(SK_'-1')")
    call disp%show( isStrAlphaNumAny(SK_'-1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAny(SK_'+1')")
    call disp%show( isStrAlphaNumAny(SK_'+1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAny(SK_'1.')")
    call disp%show( isStrAlphaNumAny(SK_'1.') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNumAny(SK_'')")
    call disp%show( isStrAlphaNumAny(SK_'') )
    call disp%skip()

end program example