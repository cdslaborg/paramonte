program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isStrAlphaAny

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isStrAlphaAny(SK_'A')")
    call disp%show( isStrAlphaAny(SK_'A') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAny(SK_'a')")
    call disp%show( isStrAlphaAny(SK_'a') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAny(SK_'aBcD')")
    call disp%show( isStrAlphaAny(SK_'aBcD') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAny(SK_'ParaMonte')")
    call disp%show( isStrAlphaAny(SK_'ParaMonte') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAny([character(4) :: 'abcd', 'ParaMonte', 'DCBA'])")
    call disp%show( isStrAlphaAny([character(4) :: 'abcd', 'ParaMonte', 'DCBA']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAny([character(9) :: 'abcd', 'ParaMonte', 'DCBA'])")
    call disp%show( isStrAlphaAny([character(9) :: 'abcd', 'ParaMonte', 'DCBA']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAny([character(4) :: ' ', '1aBc2', '#@!?'])")
    call disp%show( isStrAlphaAny([character(4) :: ' ', '1aBc2', '#@!?']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAny(SK_'-a')")
    call disp%show( isStrAlphaAny(SK_'-a') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAny(SK_'_A')")
    call disp%show( isStrAlphaAny(SK_'_A') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAny(SK_'_')")
    call disp%show( isStrAlphaAny(SK_'_') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAny(SK_' ')")
    call disp%show( isStrAlphaAny(SK_' ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAny(SK_'')")
    call disp%show( isStrAlphaAny(SK_'') )
    call disp%skip()

end program example