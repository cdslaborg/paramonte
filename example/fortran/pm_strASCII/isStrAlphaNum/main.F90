program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isStrAlphaNum

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isStrAlphaNum('1')")
    call disp%show( isStrAlphaNum('1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNum('0123')")
    call disp%show( isStrAlphaNum('0123') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNum('a01C3F')")
    call disp%show( isStrAlphaNum('a01C3F') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNum('ParaMonte')")
    call disp%show( isStrAlphaNum('ParaMonte') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNum('-1')")
    call disp%show( isStrAlphaNum('-1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNum('+1')")
    call disp%show( isStrAlphaNum('+1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNum('1.')")
    call disp%show( isStrAlphaNum('1.') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaNum('')")
    call disp%show( isStrAlphaNum('') )
    call disp%skip()

end program example