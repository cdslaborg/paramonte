program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isStrDigitAny

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isStrDigitAny(SK_'1')")
    call disp%show( isStrDigitAny(SK_'1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrDigitAny(SK_'-1')")
    call disp%show( isStrDigitAny(SK_'-1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrDigitAny(SK_'+1')")
    call disp%show( isStrDigitAny(SK_'+1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrDigitAny(SK_'1.')")
    call disp%show( isStrDigitAny(SK_'1.') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrDigitAny(SK_'')")
    call disp%show( isStrDigitAny(SK_'') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrDigitAny([character(4,SK) :: '1', 'p123', '0123'])")
    call disp%show( isStrDigitAny([character(4,SK) :: '1', 'p123', '0123']) )
    call disp%skip()

end program example