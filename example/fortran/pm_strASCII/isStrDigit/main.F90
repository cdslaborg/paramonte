program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isStrDigit

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isStrDigit(SK_'1')")
    call disp%show( isStrDigit(SK_'1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrDigit(SK_'-1')")
    call disp%show( isStrDigit(SK_'-1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrDigit(SK_'+1')")
    call disp%show( isStrDigit(SK_'+1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrDigit(SK_'1.')")
    call disp%show( isStrDigit(SK_'1.') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrDigit(SK_'')")
    call disp%show( isStrDigit(SK_'') )
    call disp%show("size(isStrDigit(SK_''))")
    call disp%show( size(isStrDigit(SK_'')) )
    call disp%skip()

end program example